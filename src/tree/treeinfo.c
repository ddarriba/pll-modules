/*
    Copyright (C) 2016 Alexey Kozlov

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Alexey Kozlov <Alexey.Kozlov@h-its.org>,
    Heidelberg Institute for Theoretical Studies,
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "pll_tree.h"

#include "../pllmod_common.h"

static int cb_set_treeinfo_pointer(pll_utree_t * node, void *data)
{
  pllmod_treeinfo_t * treeinfo = (pllmod_treeinfo_t *) data;

  ++treeinfo->counter;

  node->data = treeinfo;

  if (node->next)
  {
    node->next->data = treeinfo;
    node->next->next->data = treeinfo;
  }

  return PLL_SUCCESS;
}

/* a callback function for performing a full traversal */
static int cb_full_traversal(pll_utree_t * node)
{
  return PLL_SUCCESS;
}

/* a callback function for performing a partial traversal on invalid CLVs */
static int cb_partial_traversal(pll_utree_t * node)
{
  /* do not include tips */
  if (!node->next) return PLL_FAILURE;

  pllmod_treeinfo_t * treeinfo = (pllmod_treeinfo_t *) node->data;

  /* if clv is invalid, traverse the subtree to compute it */
  if (treeinfo->active_partition == PLLMOD_TREEINFO_PARTITION_ALL)
  {
    /* check if at least one per-partition CLV is invalid */
    size_t p;
    for (p = 0; p < treeinfo->partition_count; ++p)
      if (treeinfo->clv_valid[p][node->node_index] == 0)
        return PLL_SUCCESS;

    /* CLVs for all partitions are valid -> skip subtree */
    return PLL_FAILURE;
  }
  else
    return (treeinfo->clv_valid[treeinfo->active_partition][node->node_index] == 0);
}

static int treeinfo_partition_active(pllmod_treeinfo_t * treeinfo,
                                     int partition_index)
{
  return (treeinfo->active_partition == PLLMOD_TREEINFO_PARTITION_ALL ||
            treeinfo->active_partition == (int) partition_index);
}

PLL_EXPORT pllmod_treeinfo_t * pllmod_treeinfo_create(pll_utree_t * root,
                                                      unsigned int tips,
                                                      unsigned int partitions,
                                                      int linked_branch_lengths)
{
  /* create treeinfo instance */
  pllmod_treeinfo_t * treeinfo;

  if (!(treeinfo = (pllmod_treeinfo_t *) calloc(1, sizeof(pllmod_treeinfo_t))))
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for treeinfo\n");
    return NULL;
  }

  /* compute some derived dimensions */
  size_t inner_nodes_count = tips - 2;
  size_t nodes_count       = inner_nodes_count + tips;
  size_t branch_count      = nodes_count - 1;
  size_t pmatrix_count     = branch_count;
  size_t utree_count       = inner_nodes_count * 3 + tips;

  /* save tree pointer */
  treeinfo->root = root;

  /* traverse the tree and set treeinfo pointer into "data" field of every node */
  treeinfo->counter = 0;
  int retval = pllmod_utree_traverse_apply(root,
                                           NULL,                     /* pre */
                                           cb_set_treeinfo_pointer,  /* post */
                                           (void *) treeinfo);
  if (treeinfo->counter != nodes_count)
  {
    pllmod_set_error(PLLMOD_TREE_ERROR_INVALID_TREE_SIZE,
                     "Invalid tree size. Got %d instead of %d\n",
                     treeinfo->counter, nodes_count);
    return NULL;
  }
  if (retval != PLL_SUCCESS)
    return NULL;

  /* save dimensions & options */
  treeinfo->tip_count = tips;
  treeinfo->partition_count = partitions;
  treeinfo->linked_branches = linked_branch_lengths;

  /* allocate a buffer for storing pointers to nodes of the tree in postorder
     traversal */
  treeinfo->travbuffer = (pll_utree_t **) malloc(
                            nodes_count * sizeof(pll_utree_t *));

  /* allocate a buffer for matrix indices */
  treeinfo->matrix_indices = (unsigned int *)
                                malloc(branch_count * sizeof(unsigned int));

  /* allocate a buffer for operations (parent/child clv indices) */
  treeinfo->operations = (pll_operation_t *)
                            malloc(inner_nodes_count * sizeof(pll_operation_t));

  /* check memory allocation */
  if (!treeinfo->travbuffer || !treeinfo->matrix_indices || !treeinfo->operations)
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for treeinfo structures\n");
    return NULL;
  }

  /* allocate arrays for storing per-partition info */
  treeinfo->partitions = (pll_partition_t **) calloc(partitions, sizeof(pll_partition_t *));
  treeinfo->alphas = (double *) calloc(partitions, sizeof(double));
  treeinfo->param_indices = (unsigned int **) calloc(partitions, sizeof(unsigned int*));
  treeinfo->subst_matrix_symmetries = (int **) calloc(partitions, sizeof(int*));
  treeinfo->branch_lengths = (double **) calloc(partitions, sizeof(double*));
  treeinfo->clv_valid = (char **) calloc(partitions, sizeof(char*));
  treeinfo->pmatrix_valid = (char **) calloc(partitions, sizeof(char*));

  /* check memory allocation */
  if (!treeinfo->partitions || !treeinfo->alphas || !treeinfo->param_indices ||
      !treeinfo->subst_matrix_symmetries || !treeinfo->branch_lengths ||
      !treeinfo->clv_valid || !treeinfo->pmatrix_valid)
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for treeinfo arrays\n");
    return NULL;
  }

  size_t p;
  for (p = 0; p < partitions; ++p)
  {
    /* allocate arrays for storing the branch length */
    if (p == 0 || !linked_branch_lengths)
      treeinfo->branch_lengths[p] =
        (double *) malloc(branch_count * sizeof(double));
    else
      treeinfo->branch_lengths[p] = treeinfo->branch_lengths[0];

    /* allocate invalidation arrays */
    treeinfo->clv_valid[p] = (char *) calloc(utree_count, sizeof(char));
    treeinfo->pmatrix_valid[p] = (char *) calloc(pmatrix_count, sizeof(char));

    /* check memory allocation */
    if (!(treeinfo->branch_lengths[p] &&
          treeinfo->clv_valid[p] &&
          treeinfo->pmatrix_valid[p]))
    {
      pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                  "Cannot allocate memory for arrays for partition %d\n",
                  p);
      return NULL;
    }
  }

  /* by default, work with all partitions */
  treeinfo->active_partition = PLLMOD_TREEINFO_PARTITION_ALL;

  return treeinfo;
}

PLL_EXPORT int pllmod_treeinfo_init_partition(pllmod_treeinfo_t * treeinfo,
                                           unsigned int partition_index,
                                           pll_partition_t * partition,
                                           double alpha,
                                           const unsigned int * param_indices,
                                           const int * subst_matrix_symmetries)
{
  if (!treeinfo)
  {
    pllmod_set_error(PLL_ERROR_PARAM_INVALID,
              "Treeinfo structure is NULL\n");
    return PLL_FAILURE;
  }
  else if (partition_index >= treeinfo->partition_count)
  {
    pllmod_set_error(PLL_ERROR_PARAM_INVALID,
              "Partition %d is out of bounds\n", partition_index);
    return PLL_FAILURE;
  }

  treeinfo->partitions[partition_index] = partition;
  treeinfo->alphas[partition_index] = alpha;

  /* allocate param_indices array and initialize it to all 0s,
   * i.e. per default, all rate categories will use
   * the same substitution matrix and same base frequencies */
  treeinfo->param_indices[partition_index] =
    (unsigned int *) calloc(partition->rate_cats, sizeof(unsigned int));

  /* check memory allocation */
  if (!treeinfo->param_indices[partition_index])
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for parameter indices\n");
    return PLL_FAILURE;
  }

  /* if param_indices were provided, use them instead of default */
  if (param_indices)
    memcpy(treeinfo->param_indices[partition_index],
           param_indices,
           partition->rate_cats * sizeof(unsigned int));

  /* copy substitution rate matrix symmetries, if any */
  if (subst_matrix_symmetries)
  {
    const size_t symm_size = (partition->states * (partition->states - 1) / 2)
                              * sizeof(int);
    treeinfo->subst_matrix_symmetries[partition_index] =
                               (int *) malloc(symm_size);

    /* check memory allocation */
    if (!treeinfo->subst_matrix_symmetries[partition_index])
    {
      pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                    "Cannot allocate memory for substitution scheme\n");
      return PLL_FAILURE;
    }

    memcpy(treeinfo->subst_matrix_symmetries[partition_index],
           subst_matrix_symmetries,
           symm_size);
  }
  else
    treeinfo->subst_matrix_symmetries[partition_index] = NULL;

  return PLL_SUCCESS;
}

PLL_EXPORT int pllmod_treeinfo_set_active_partition(pllmod_treeinfo_t * treeinfo,
                                                    int partition_index)
{
  if (partition_index != PLLMOD_TREEINFO_PARTITION_ALL &&
      partition_index >= (int) treeinfo->partition_count)
  {
    pllmod_set_error(PLL_ERROR_PARAM_INVALID,
              "Partition %d is out of bounds\n", partition_index);
    return PLL_FAILURE;
  }
  else
  {
    treeinfo->active_partition = partition_index;
    return PLL_SUCCESS;
  }
}

PLL_EXPORT void pllmod_treeinfo_set_root(pllmod_treeinfo_t * treeinfo,
                                         pll_utree_t * root)
{
  assert(treeinfo && root);

  /* root must be an inner node! */
  treeinfo->root = pllmod_utree_is_tip(root) ? root->back : root;
}

PLL_EXPORT void pllmod_treeinfo_set_branch_length(pllmod_treeinfo_t * treeinfo,
                                                  pll_utree_t * edge,
                                                  double length)
{
  pllmod_utree_set_length(edge, length);

  /* invalidate p-matrices */
  pllmod_treeinfo_invalidate_pmatrix(treeinfo, edge);

  /* invalidate CLVs */
  pllmod_treeinfo_invalidate_clv(treeinfo, edge->next);
  pllmod_treeinfo_invalidate_clv(treeinfo, edge->next->next);
  pllmod_treeinfo_invalidate_clv(treeinfo, edge->back->next);
  pllmod_treeinfo_invalidate_clv(treeinfo, edge->back->next->next);
}

PLL_EXPORT int pllmod_treeinfo_destroy_partition(pllmod_treeinfo_t * treeinfo,
                                                  unsigned int partition_index)
{
  if (!treeinfo)
  {
    pllmod_set_error(PLL_ERROR_PARAM_INVALID,
              "Treeinfo structure is NULL\n");
    return PLL_FAILURE;
  }
  else if (partition_index >= treeinfo->partition_count)
  {
    pllmod_set_error(PLL_ERROR_PARAM_INVALID,
              "Partition %d is out of bounds\n", partition_index);
    return PLL_FAILURE;
  }

  if (treeinfo->param_indices[partition_index])
  {
    free(treeinfo->param_indices[partition_index]);
    treeinfo->param_indices[partition_index] = NULL;
  }

  if (treeinfo->subst_matrix_symmetries[partition_index])
  {
    free(treeinfo->subst_matrix_symmetries[partition_index]);
    treeinfo->subst_matrix_symmetries[partition_index] = NULL;
  }

  return PLL_SUCCESS;
}

PLL_EXPORT void pllmod_treeinfo_destroy(pllmod_treeinfo_t * treeinfo)
{
  if (!treeinfo) return;

  /* deallocate traversal buffer, branch lengths array, matrix indices
     array and operations */
  free(treeinfo->travbuffer);
  free(treeinfo->matrix_indices);
  free(treeinfo->operations);

  /* destroy all structures allocated for the concrete PLL partition instance */
  size_t p;
  for (p = 0; p < treeinfo->partition_count; ++p)
  {
    free(treeinfo->clv_valid[p]);
    free(treeinfo->pmatrix_valid[p]);
    free(treeinfo->param_indices[p]);

    if (p == 0 || !treeinfo->linked_branches)
      free(treeinfo->branch_lengths[p]);
    if (treeinfo->subst_matrix_symmetries)
      free(treeinfo->subst_matrix_symmetries[p]);
  }

  if(treeinfo->subst_matrix_symmetries)
    free(treeinfo->subst_matrix_symmetries);

  /* free invalidation arrays */
  free(treeinfo->clv_valid);
  free(treeinfo->pmatrix_valid);

  /* free alpha and param_indices arrays */
  free(treeinfo->alphas);
  free(treeinfo->param_indices);
  free(treeinfo->branch_lengths);

  /* deallocate partition array */
  free(treeinfo->partitions);

  /* finally, deallocate treeinfo object itself */
  free(treeinfo);
}

PLL_EXPORT int pllmod_treeinfo_update_prob_matrices(pllmod_treeinfo_t * treeinfo,
                                                    int update_all)
{
  size_t p, m;
  size_t updated = 0;
  size_t pmatrix_count = 2 * treeinfo->tip_count - 3;

  for (p = 0; p < treeinfo->partition_count; ++p)
  {
    /* only selected partitioned will be affected */
    if (treeinfo_partition_active(treeinfo, p))
    {

      for (m = 0; m < pmatrix_count; ++m)
      {
        const unsigned int matrix_index = treeinfo->matrix_indices[m];

        if (treeinfo->pmatrix_valid[p][matrix_index] && !update_all)
          continue;

        pll_update_prob_matrices (treeinfo->partitions[p],
                                  treeinfo->param_indices[p],
                                  &matrix_index,
                                  &treeinfo->branch_lengths[p][m],
                                  1);

        treeinfo->pmatrix_valid[p][matrix_index] = 1;
        updated++;
      }
    }
  }

  return PLL_SUCCESS;
}

PLL_EXPORT void pllmod_treeinfo_invalidate_all(pllmod_treeinfo_t * treeinfo)
{
  size_t p, m;
  size_t clv_count = treeinfo->tip_count + (treeinfo->tip_count - 2) * 3;
  size_t pmatrix_count = 2 * treeinfo->tip_count - 3;

  for (p = 0; p < treeinfo->partition_count; ++p)
  {
    /* only selected partitioned will be affected */
    if (treeinfo_partition_active(treeinfo, p))
    {
      for (m = 0; m < pmatrix_count; ++m)
        treeinfo->pmatrix_valid[p][m] = 0;

      for (m = 0; m < clv_count; ++m)
        treeinfo->clv_valid[p][m] = 0;
    }
  }
}

PLL_EXPORT int pllmod_treeinfo_validate_clvs(pllmod_treeinfo_t * treeinfo,
                                             pll_utree_t ** travbuffer,
                                             unsigned int travbuffer_size)
{
  size_t p;
  for (p = 0; p < treeinfo->partition_count; ++p)
  {
    /* only selected partitioned will be affected */
    if (treeinfo_partition_active(treeinfo, p))
    {
      unsigned int i;
      for (i = 0; i < travbuffer_size; ++i)
      {
        const pll_utree_t * node = travbuffer[i];
        if (node->next)
        {
          treeinfo->clv_valid[p][node->node_index] = 1;

          /* since we have only 1 CLV vector per inner node,
           * we must invalidate CLVs for other 2 directions */
          treeinfo->clv_valid[p][node->next->node_index] = 0;
          treeinfo->clv_valid[p][node->next->next->node_index] = 0;
        }
      }
    }
  }

  return PLL_SUCCESS;
}

PLL_EXPORT void pllmod_treeinfo_invalidate_pmatrix(pllmod_treeinfo_t * treeinfo,
                                                   const pll_utree_t * edge)
{
  size_t p;
  for (p = 0; p < treeinfo->partition_count; ++p)
  {
    if (treeinfo_partition_active(treeinfo, p))
      treeinfo->pmatrix_valid[p][edge->pmatrix_index] = 0;
  }
}

PLL_EXPORT void pllmod_treeinfo_invalidate_clv(pllmod_treeinfo_t * treeinfo,
                                               const pll_utree_t * edge)
{
  size_t p;
  for (p = 0; p < treeinfo->partition_count; ++p)
  {
    if (treeinfo_partition_active(treeinfo, p))
      treeinfo->clv_valid[p][edge->node_index] = 0;
  }
}

PLL_EXPORT double pllmod_treeinfo_compute_loglh(pllmod_treeinfo_t * treeinfo,
                                                int incremental)
{
  /* tree root must be an inner node! */
  assert(!pllmod_utree_is_tip(treeinfo->root));

  const double LOGLH_NONE = NAN;

  double total_loglh = 0.0;

  const int old_active_partition = treeinfo->active_partition;

  /* iterate over all partitions (we assume that traversal is the same) */
  size_t p;
  for (p = 0; p < treeinfo->partition_count; ++p)
  {
    /* all subsequent operation will affect current partition only */
    pllmod_treeinfo_set_active_partition(treeinfo, p);

    /* perform a FULL postorder traversal of the unrooted tree */
    unsigned int traversal_size;
    if (!pll_utree_traverse(treeinfo->root,
                            cb_full_traversal,
                            treeinfo->travbuffer,
                            &traversal_size))
      return LOGLH_NONE;


    /* given the computed traversal descriptor, generate the operations
       structure, and the corresponding probability matrix indices that
       may need recomputing */
    unsigned int matrix_count, ops_count;
    pll_utree_create_operations(treeinfo->travbuffer,
                                traversal_size,
                                treeinfo->branch_lengths[p],
                                treeinfo->matrix_indices,
                                treeinfo->operations,
                                &matrix_count,
                                &ops_count);

    pllmod_treeinfo_update_prob_matrices(treeinfo, !incremental);

    if (incremental)
    {
      /* compute partial traversal and update only invalid CLVs */
      if (!pll_utree_traverse(treeinfo->root,
                              cb_partial_traversal,
                              treeinfo->travbuffer,
                              &traversal_size))
        return LOGLH_NONE;

      /* create operations based on partial traversal */
      pll_utree_create_operations(treeinfo->travbuffer,
                                  traversal_size,
                                  treeinfo->branch_lengths[p],
                                  treeinfo->matrix_indices,
                                  treeinfo->operations,
                                  &matrix_count,
                                  &ops_count);
    }

    treeinfo->counter += ops_count;

    /* use the operations array to compute all ops_count inner CLVs. Operations
       will be carried out sequentially starting from operation 0 towards
       ops_count-1 */
    pll_update_partials(treeinfo->partitions[p],
                        treeinfo->operations,
                        ops_count);

    pllmod_treeinfo_validate_clvs(treeinfo,
                                  treeinfo->travbuffer,
                                  traversal_size);

    /* compute the likelihood on an edge of the unrooted tree by specifying
       the CLV indices at the two end-point of the branch, the probability
       matrix index for the concrete branch length, and the index of the model
       of whose frequency vector is to be used */
    const double part_loglh = pll_compute_edge_loglikelihood(
                                            treeinfo->partitions[p],
                                            treeinfo->root->clv_index,
                                            treeinfo->root->scaler_index,
                                            treeinfo->root->back->clv_index,
                                            treeinfo->root->back->scaler_index,
                                            treeinfo->root->pmatrix_index,
                                            treeinfo->param_indices[p],
                                            NULL);

    /* accumalate loglh by summing up over all the partitions */
    total_loglh += part_loglh;
  }

  /* restore original active partition */
  pllmod_treeinfo_set_active_partition(treeinfo, old_active_partition);

  assert(total_loglh < 0.);

  return total_loglh;
}
