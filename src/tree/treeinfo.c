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

static int treeinfo_check_tree(pllmod_treeinfo_t * treeinfo,
                               pll_utree_t * tree);
static int treeinfo_init_tree(pllmod_treeinfo_t * treeinfo);

/* a callback function for performing a full traversal */
static int cb_full_traversal(pll_unode_t * node)
{
  PLLMOD_UNUSED(node);
  return PLL_SUCCESS;
}

/* a callback function for performing a partial traversal on invalid CLVs */
static int cb_partial_traversal(pll_unode_t * node)
{
  /* do not include tips */
  if (!node->next) return PLL_FAILURE;

  pllmod_treeinfo_t * treeinfo = (pllmod_treeinfo_t *) node->data;

  /* if clv is invalid, traverse the subtree to compute it */
  if (treeinfo->active_partition == PLLMOD_TREEINFO_PARTITION_ALL)
  {
    /* check if at least one per-partition CLV is invalid */
    for (unsigned int i = 0; i < treeinfo->init_partition_count; ++i)
    {
      unsigned int p = treeinfo->init_partition_idx[i];
      if (treeinfo->clv_valid[p][node->node_index] == 0)
        return PLL_SUCCESS;
    }

    /* CLVs for all partitions are valid -> skip subtree */
    return PLL_FAILURE;
  }
  else
    return (treeinfo->clv_valid[treeinfo->active_partition][node->node_index] == 0);
}

static int treeinfo_partition_active(pllmod_treeinfo_t * treeinfo,
                                     unsigned int partition_index)
{
  return (treeinfo->active_partition == PLLMOD_TREEINFO_PARTITION_ALL ||
          treeinfo->active_partition == (int) partition_index);
}

PLL_EXPORT pllmod_treeinfo_t * pllmod_treeinfo_create(pll_unode_t * root,
                                                      unsigned int tips,
                                                      unsigned int partitions,
                                                      int brlen_linkage)
{
  /* create treeinfo instance */
  pllmod_treeinfo_t * treeinfo;

  if (!(treeinfo = (pllmod_treeinfo_t *) calloc(1, sizeof(pllmod_treeinfo_t))))
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for treeinfo\n");
    return NULL;
  }

  /* save dimensions & options */
  treeinfo->tip_count = tips;
  treeinfo->partition_count = partitions;
  treeinfo->brlen_linkage = brlen_linkage;

  /* create pll_utree structure and store it in treeinfo */
  treeinfo->tree = pll_utree_wraptree(root, tips);

  if (!treeinfo->tree || !treeinfo_check_tree(treeinfo, treeinfo->tree))
  {
    assert(pll_errno);
    free(treeinfo);
    return NULL;
  }

  /* compute some derived dimensions */
  unsigned int inner_nodes_count = treeinfo->tree->inner_count;
  unsigned int nodes_count       = inner_nodes_count + tips;
  unsigned int branch_count      = treeinfo->tree->edge_count;
  treeinfo->subnode_count        = tips + 3 * inner_nodes_count;

  /* allocate a buffer for storing pointers to nodes of the tree in postorder
     traversal */
  treeinfo->travbuffer = (pll_unode_t **) malloc(
                            nodes_count * sizeof(pll_unode_t *));

  /* allocate a buffer for matrix indices */
  treeinfo->matrix_indices = (unsigned int *)
                                malloc(branch_count * sizeof(unsigned int));

  /* allocate a buffer for operations (parent/child clv indices) */
  treeinfo->operations = (pll_operation_t *)
                            malloc(inner_nodes_count * sizeof(pll_operation_t));

  treeinfo->subnodes = (pll_unode_t **) malloc(
                          treeinfo->subnode_count * sizeof(pll_unode_t *));

  /* check memory allocation */
  if (!treeinfo->travbuffer || !treeinfo->matrix_indices ||
      !treeinfo->operations || !treeinfo->subnodes)
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for treeinfo structures\n");
    return NULL;
  }

  /* allocate arrays for storing per-partition info */
  treeinfo->partitions = (pll_partition_t **) calloc(partitions, sizeof(pll_partition_t *));
  treeinfo->params_to_optimize = (int *) calloc(partitions, sizeof(int));
  treeinfo->alphas = (double *) calloc(partitions, sizeof(double));
  treeinfo->gamma_mode = (int *) calloc(partitions, sizeof(int));
  treeinfo->param_indices = (unsigned int **) calloc(partitions, sizeof(unsigned int*));
  treeinfo->subst_matrix_symmetries = (int **) calloc(partitions, sizeof(int*));
  treeinfo->branch_lengths = (double **) calloc(partitions, sizeof(double*));
  treeinfo->deriv_precomp = (double **) calloc(partitions, sizeof(double*));
  treeinfo->clv_valid = (char **) calloc(partitions, sizeof(char*));
  treeinfo->pmatrix_valid = (char **) calloc(partitions, sizeof(char*));
  treeinfo->partition_loglh = (double *) calloc(partitions, sizeof(double));

  treeinfo->init_partition_count = 0;
  treeinfo->init_partition_idx = (unsigned int *) calloc(partitions, sizeof(unsigned int));
  treeinfo->init_partitions = (pll_partition_t **) calloc(partitions, sizeof(pll_partition_t *));

  /* allocate array for storing linked/average branch lengths */
  treeinfo->linked_branch_lengths = (double *) malloc(branch_count * sizeof(double));

  /* allocate branch length scalers if needed */
  if (brlen_linkage == PLLMOD_COMMON_BRLEN_SCALED)
    treeinfo->brlen_scalers = (double *) calloc(partitions, sizeof(double));
  else
    treeinfo->brlen_scalers = NULL;

  /* check memory allocation */
  if (!treeinfo->partitions || !treeinfo->alphas || !treeinfo->param_indices ||
      !treeinfo->subst_matrix_symmetries || !treeinfo->branch_lengths ||
      !treeinfo->deriv_precomp || !treeinfo->clv_valid || !treeinfo->pmatrix_valid ||
      !treeinfo->linked_branch_lengths || !treeinfo->partition_loglh ||
      !treeinfo->gamma_mode || !treeinfo->init_partition_idx ||
      (brlen_linkage == PLLMOD_COMMON_BRLEN_SCALED && !treeinfo->brlen_scalers))
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for treeinfo arrays\n");
    return NULL;
  }

  unsigned int p;
  for (p = 0; p < partitions; ++p)
  {
    /* use mean GAMMA rates per default */
    treeinfo->gamma_mode[p] = PLL_GAMMA_RATES_MEAN;

    /* allocate arrays for storing the per-partition branch lengths */
    if (brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED)
    {
      treeinfo->branch_lengths[p] =
          (double *) malloc(branch_count * sizeof(double));
    }
    else
      treeinfo->branch_lengths[p] = treeinfo->linked_branch_lengths;

    /* initialize all branch length scalers to 1 */
    if (treeinfo->brlen_scalers)
      treeinfo->brlen_scalers[p] = 1.;

    /* check memory allocation */
    if (!treeinfo->branch_lengths[p])
    {
      pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                  "Cannot allocate memory for arrays for partition %d\n",
                  p);
      return NULL;
    }
  }

  /* by default, work with all partitions */
  treeinfo->active_partition = PLLMOD_TREEINFO_PARTITION_ALL;

  /* needs to be here since we use some of the arrays allocated above */
  if (!treeinfo_init_tree(treeinfo))
  {
    pllmod_treeinfo_destroy(treeinfo);
    assert(pll_errno);
    return NULL;
  }

  assert(treeinfo->tree && treeinfo->tree->tip_count == tips);

  return treeinfo;
}

PLL_EXPORT
int pllmod_treeinfo_set_parallel_context(pllmod_treeinfo_t * treeinfo,
                                         void * parallel_context,
                                         void (*parallel_reduce_cb)(void *,
                                                                    double *,
                                                                    size_t,
                                                                    int))
{
  treeinfo->parallel_context = parallel_context;
  treeinfo->parallel_reduce_cb = parallel_reduce_cb;

  return PLL_SUCCESS;
}


PLL_EXPORT int pllmod_treeinfo_init_partition(pllmod_treeinfo_t * treeinfo,
                                           unsigned int partition_index,
                                           pll_partition_t * partition,
                                           int params_to_optimize,
                                           int gamma_mode,
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
  else if (treeinfo->partitions[partition_index])
  {
    pllmod_set_error(PLL_ERROR_PARAM_INVALID,
              "Partition %d is already initialized\n", partition_index);
    return PLL_FAILURE;
  }

  unsigned int local_partition_index = treeinfo->init_partition_count++;
  treeinfo->partitions[partition_index] = partition;
  treeinfo->init_partitions[local_partition_index] = partition;
  treeinfo->init_partition_idx[local_partition_index] = partition_index;
  treeinfo->params_to_optimize[partition_index] = params_to_optimize;
  treeinfo->gamma_mode[partition_index] = gamma_mode;
  treeinfo->alphas[partition_index] = alpha;

  /* compute some derived dimensions */
  unsigned int inner_nodes_count = treeinfo->tip_count - 2;
  unsigned int nodes_count       = inner_nodes_count + treeinfo->tip_count;
  unsigned int branch_count      = nodes_count - 1;
  unsigned int pmatrix_count     = branch_count;
  unsigned int utree_count       = inner_nodes_count * 3 + treeinfo->tip_count;

  /* allocate invalidation arrays */
  treeinfo->clv_valid[partition_index] =
      (char *) calloc(utree_count, sizeof(char));
  treeinfo->pmatrix_valid[partition_index] = (
      char *) calloc(pmatrix_count, sizeof(char));

  /* check memory allocation */
  if (!treeinfo->clv_valid[partition_index] ||
      !treeinfo->pmatrix_valid[partition_index])
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for parameter indices\n");
    return PLL_FAILURE;
  }

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
    const unsigned int symm_size =
                            (partition->states * (partition->states - 1) / 2)
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

  /* allocate memory for derivative precomputation table */
  unsigned int sites_alloc = partition->sites;
  if (partition->attributes & PLL_ATTRIB_AB_FLAG)
    sites_alloc += partition->states;
  unsigned int precomp_size =
      sites_alloc * partition->rate_cats * partition->states_padded;

  treeinfo->deriv_precomp[partition_index] = (double *) pll_aligned_alloc(
      precomp_size * sizeof(double), partition->alignment);

  if (!treeinfo->deriv_precomp[partition_index])
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                  "Cannot allocate memory for derivative buffers\n");
    return PLL_FAILURE;
  }

  memset(treeinfo->deriv_precomp[partition_index], 0, precomp_size * sizeof(double));

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

PLL_EXPORT int pllmod_treeinfo_set_root(pllmod_treeinfo_t * treeinfo,
                                         pll_unode_t * root)
{
  if (!treeinfo || !root || root->data != (void *) treeinfo)
  {
    pllmod_set_error(PLL_ERROR_PARAM_INVALID, "Invalid root node!\n");
    return PLL_FAILURE;
  }

  /* root must be an inner node! */
  treeinfo->root = pllmod_utree_is_tip(root) ? root->back : root;
  treeinfo->tree->vroot = treeinfo->root;

  return PLL_SUCCESS;
}

PLL_EXPORT
int pllmod_treeinfo_get_branch_length_all(const pllmod_treeinfo_t * treeinfo,
                                          const pll_unode_t * edge,
                                          double * lengths)
{
  unsigned int pmatrix_index = edge->pmatrix_index;

  if (treeinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED)
  {
    for (unsigned int p = 0; p < treeinfo->partition_count; ++p)
    {
      if (treeinfo->partitions[p])
        *lengths++ = treeinfo->branch_lengths[p][pmatrix_index];
    }
  }
  else
    lengths[0] = treeinfo->branch_lengths[0][pmatrix_index];

  return PLL_SUCCESS;
}

PLL_EXPORT int pllmod_treeinfo_set_branch_length(pllmod_treeinfo_t * treeinfo,
                                                  pll_unode_t * edge,
                                                  double length)
{
  assert (treeinfo->brlen_linkage != PLLMOD_COMMON_BRLEN_UNLINKED);
  return pllmod_treeinfo_set_branch_length_all(treeinfo, edge, &length);
}

PLL_EXPORT
int pllmod_treeinfo_set_branch_length_all(pllmod_treeinfo_t * treeinfo,
                                          pll_unode_t * edge,
                                          const double * lengths)
{
  unsigned int pmatrix_index = edge->pmatrix_index;
  if (treeinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED)
  {
    for (unsigned int p = 0; p < treeinfo->partition_count; ++p)
    {
      if (treeinfo->partitions[p])
        treeinfo->branch_lengths[p][pmatrix_index] = *lengths++;
    }
  }
  else
  {
    treeinfo->branch_lengths[0][pmatrix_index] = lengths[0];
    pllmod_utree_set_length(edge, lengths[0]);
  }

#if 0
  pllmod_treeinfo_set_active_partition(treeinfo, PLLMOD_TREEINFO_PARTITION_ALL);

  /* invalidate p-matrices */
  pllmod_treeinfo_invalidate_pmatrix(treeinfo, edge);

  /* invalidate CLVs */
  if (edge->next)
  {
    pllmod_treeinfo_invalidate_clv(treeinfo, edge->next);
    pllmod_treeinfo_invalidate_clv(treeinfo, edge->next->next);
  }
  if (edge->back->next)
  {
    pllmod_treeinfo_invalidate_clv(treeinfo, edge->back->next);
    pllmod_treeinfo_invalidate_clv(treeinfo, edge->back->next->next);
  }
#endif

  return PLL_SUCCESS;
}

PLL_EXPORT
int pllmod_treeinfo_set_branch_length_partition(pllmod_treeinfo_t * treeinfo,
                                                pll_unode_t * edge,
                                                int partition_index,
                                                double length)
{
  unsigned int pmatrix_index = edge->pmatrix_index;
  const int old_active_partition = treeinfo->active_partition;

  if (!pllmod_treeinfo_set_active_partition(treeinfo, partition_index))
      return PLL_FAILURE;

  if (partition_index != PLLMOD_TREEINFO_PARTITION_ALL)
    treeinfo->branch_lengths[partition_index][pmatrix_index] = length;
  else if (treeinfo->brlen_linkage != PLLMOD_COMMON_BRLEN_UNLINKED)
  {
    pllmod_utree_set_length(edge, length);
    treeinfo->branch_lengths[0][pmatrix_index] = length;
  }
  else
  {
    for (unsigned int p = 0; p < treeinfo->partition_count; ++p)
    {
      if (treeinfo->partitions[p])
        treeinfo->branch_lengths[p][pmatrix_index] = length;
    }
  }

  treeinfo->active_partition = old_active_partition;

#if 0
  /* invalidate p-matrices */
  pllmod_treeinfo_invalidate_pmatrix(treeinfo, edge);

  /* invalidate CLVs */
  if (edge->next)
  {
    pllmod_treeinfo_invalidate_clv(treeinfo, edge->next);
    pllmod_treeinfo_invalidate_clv(treeinfo, edge->next->next);
  }
  if (edge->back->next)
  {
    pllmod_treeinfo_invalidate_clv(treeinfo, edge->back->next);
    pllmod_treeinfo_invalidate_clv(treeinfo, edge->back->next->next);
  }
#endif

  return PLL_SUCCESS;
}

PLL_EXPORT
pll_utree_t * pllmod_treeinfo_get_partition_tree(const pllmod_treeinfo_t * treeinfo,
                                                 int partition_index)
{
  pll_utree_t * ptree = pll_utree_clone(treeinfo->tree);

  // set partition-specific branch lengths
  if (treeinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED)
  {
    unsigned int node_count = ptree->tip_count + ptree->inner_count;
    unsigned int edge_count = ptree->edge_count;
    const double * brlens = treeinfo->branch_lengths[partition_index];
    for (unsigned int i = 0; i < node_count; ++i)
    {
      pll_unode_t * snode = ptree->nodes[i];
      do
      {
        unsigned int pmat_idx = snode->pmatrix_index;
        if (pmat_idx > edge_count)
        {
          pllmod_set_error(PLLMOD_TREE_ERROR_INVALID_TREE,
                           "p-matrix index out of bounds (%u). "
                           "treeinfo structure requires that each branch "
                           "is assigned a unique p-matrix index "
                           "between 0 and branch_count-1\n",
                           pmat_idx);
          return NULL;
        }
        snode->length = brlens[pmat_idx];
        snode = snode->next;
      }
      while (snode && snode != ptree->nodes[i]);
    }
  }

  return ptree;
}

PLL_EXPORT
pllmod_treeinfo_topology_t * pllmod_treeinfo_get_topology(const pllmod_treeinfo_t * treeinfo,
                                                          pllmod_treeinfo_topology_t * topol)
{
  unsigned int brlen_set_count =
      (treeinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED) ?
          treeinfo->init_partition_count : 0;

  if (!topol)
  {
    topol =
        (pllmod_treeinfo_topology_t *) calloc(1, sizeof(pllmod_treeinfo_topology_t));
    if (!topol)
    {
      pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                    "Cannot allocate memory for topology structure\n");
      return PLL_FAILURE;
    }

    topol->edge_count = treeinfo->tree->edge_count;
    topol->brlen_set_count = brlen_set_count;
    topol->root_index = treeinfo->root->node_index;
    topol->edges = (pllmod_treeinfo_edge_t *) calloc(topol->edge_count,
                                                     sizeof(pllmod_treeinfo_edge_t));
    if (!topol->edges)
    {
      pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                    "Cannot allocate memory for topology buffers\n");
      pllmod_treeinfo_destroy_topology(topol);
      return PLL_FAILURE;
    }

    if (topol->brlen_set_count)
    {
      topol->branch_lengths = (double **) calloc(topol->brlen_set_count,
                                                 sizeof(double *));
      if (!topol->branch_lengths)
      {
        pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                      "Cannot allocate memory for branch length buffers\n");
        pllmod_treeinfo_destroy_topology(topol);
        return PLL_FAILURE;
      }

      for (unsigned int i = 0; i < topol->brlen_set_count; ++i)
      {
        topol->branch_lengths[i] = (double *) calloc(topol->edge_count,
                                                     sizeof(double));
        if (!topol->branch_lengths[i])
        {
          pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                        "Cannot allocate memory for branch length buffers\n");
          pllmod_treeinfo_destroy_topology(topol);
          return PLL_FAILURE;
        }
      }
    }
  }
  else if (treeinfo->tree->edge_count != topol->edge_count ||
           brlen_set_count != topol->brlen_set_count)
  {
    pllmod_set_error(PLL_ERROR_PARAM_INVALID,
                     "Incompatible treeinfo_topology structure!\n");
    return PLL_FAILURE;
  }

  // save topology as a list of edges
  unsigned int edge_num = 0;
  for (unsigned int i = 0; i < treeinfo->subnode_count; ++i)
  {
    const pll_unode_t * snode = treeinfo->subnodes[i];
    if (snode->node_index < snode->back->node_index)
    {
      topol->edges[edge_num].pmatrix_index = snode->pmatrix_index;
      topol->edges[edge_num].left_index = snode->node_index;
      topol->edges[edge_num].right_index = snode->back->node_index;
      topol->edges[edge_num].brlen = snode->length;
      edge_num++;
    }
  }
  assert(edge_num == topol->edge_count);

  // save unlinked branch lengths
  for (unsigned int i = 0; i < topol->brlen_set_count; ++i)
  {
    unsigned int p = treeinfo->init_partition_idx[i];
    memcpy(topol->branch_lengths[i], treeinfo->branch_lengths[p],
            topol->edge_count * sizeof(double));
  }

  return topol;
}

PLL_EXPORT
int pllmod_treeinfo_set_topology(pllmod_treeinfo_t * treeinfo,
                                 const pllmod_treeinfo_topology_t * topol)
{
  unsigned int brlen_set_count;

  if (!treeinfo || !topol)
  {
    pllmod_set_error(PLL_ERROR_PARAM_INVALID, "treeinfo and/or topology is empty!\n");
    return PLL_FAILURE;
  }
  if (treeinfo->tree->edge_count != topol->edge_count)
  {
    pllmod_set_error(PLL_ERROR_PARAM_INVALID,
                     "Incompatible topology: edge count differs!\n");
    return PLL_FAILURE;
  }

  brlen_set_count = (treeinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED) ?
                     treeinfo->init_partition_count : 0;
  if (brlen_set_count != topol->brlen_set_count)
  {
    pllmod_set_error(PLL_ERROR_PARAM_INVALID,
                     "Incompatible topology: brlen set count differs!!\n");
    return PLL_FAILURE;
  }

  // re-connect branches and reset pmatrix indices
  for (unsigned int i = 0; i < topol->edge_count; ++i)
  {
    const pllmod_treeinfo_edge_t * edge = &topol->edges[i];
    pll_unode_t * left_node = treeinfo->subnodes[edge->left_index];
    pll_unode_t * right_node = treeinfo->subnodes[edge->right_index];
    pllmod_utree_connect_nodes(left_node, right_node, edge->brlen);
    left_node->pmatrix_index = right_node->pmatrix_index = edge->pmatrix_index;

//    printf("load edge: %u %u %u %f\n", left_node->node_index, right_node->node_index,
//           left_node->pmatrix_index, left_node->length);

    // linked/scaled brlens -> set a global value for all partitions
    if (!topol->branch_lengths)
      treeinfo->branch_lengths[0][edge->pmatrix_index] = edge->brlen;
  }

  // restore unlinked branch lengths
  if (topol->branch_lengths)
  {
    for (unsigned int i = 0; i < topol->brlen_set_count; ++i)
    {
      unsigned int p = treeinfo->init_partition_idx[i];
      memcpy(treeinfo->branch_lengths[p], topol->branch_lengths[i],
             topol->edge_count * sizeof(double));
    }
  }

  // restore virtual root
  treeinfo->root = treeinfo->tree->vroot = treeinfo->subnodes[topol->root_index];
  assert(treeinfo->root->node_index == topol->root_index);

  pllmod_treeinfo_invalidate_all(treeinfo);

  return PLL_SUCCESS;
}

PLL_EXPORT
int pllmod_treeinfo_destroy_topology(pllmod_treeinfo_topology_t * topol)
{
  if (topol)
  {
    if (topol->edges)
      free(topol->edges);
    if (topol->branch_lengths)
    {
      for (unsigned int i = 0; i < topol->brlen_set_count; ++i)
        free(topol->branch_lengths[i]);
      free(topol->branch_lengths);
    }
    free(topol);
  }
  return PLL_SUCCESS;
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

  if (treeinfo->clv_valid[partition_index])
  {
    free(treeinfo->clv_valid[partition_index]);
    treeinfo->clv_valid[partition_index] = NULL;
  }

  if (treeinfo->pmatrix_valid[partition_index])
  {
    free(treeinfo->pmatrix_valid[partition_index]);
    treeinfo->pmatrix_valid[partition_index] = NULL;
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

  if (treeinfo->deriv_precomp[partition_index])
  {
    free(treeinfo->deriv_precomp[partition_index]);
    treeinfo->deriv_precomp[partition_index] = NULL;
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
  free(treeinfo->subnodes);

  /* destroy all structures allocated for the concrete PLL partition instance */
  unsigned int p;
  for (p = 0; p < treeinfo->partition_count; ++p)
  {
    if (treeinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED)
      free(treeinfo->branch_lengths[p]);

    pllmod_treeinfo_destroy_partition(treeinfo, p);
  }

  if(treeinfo->subst_matrix_symmetries)
    free(treeinfo->subst_matrix_symmetries);

  if(treeinfo->constraint)
    free(treeinfo->constraint);

  /* free invalidation arrays */
  free(treeinfo->clv_valid);
  free(treeinfo->pmatrix_valid);

  free(treeinfo->linked_branch_lengths);

  /* free alpha and param_indices arrays */
  free(treeinfo->params_to_optimize);
  free(treeinfo->alphas);
  free(treeinfo->gamma_mode);
  free(treeinfo->param_indices);
  free(treeinfo->branch_lengths);
  free(treeinfo->partition_loglh);
  free(treeinfo->deriv_precomp);

  if(treeinfo->brlen_scalers)
    free(treeinfo->brlen_scalers);

  /* deallocate partition array */
  free(treeinfo->partitions);
  free(treeinfo->init_partitions);
  free(treeinfo->init_partition_idx);

  if (treeinfo->tree)
  {
    free(treeinfo->tree->nodes);
    free(treeinfo->tree);
  }

  /* finally, deallocate treeinfo object itself */
  free(treeinfo);
}

PLL_EXPORT int pllmod_treeinfo_update_prob_matrices(pllmod_treeinfo_t * treeinfo,
                                                    int update_all)
{
  unsigned int p, m;
  unsigned int updated = 0;
  unsigned int pmatrix_count = treeinfo->tree->edge_count;

  for (p = 0; p < treeinfo->partition_count; ++p)
  {
    /* only selected partitioned will be affected */
    if (treeinfo_partition_active(treeinfo, p) && treeinfo->partitions[p])
    {

      for (m = 0; m < pmatrix_count; ++m)
      {
        if (treeinfo->pmatrix_valid[p][m] && !update_all)
          continue;

        double p_brlen = treeinfo->branch_lengths[p][m];
        if (treeinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_SCALED)
          p_brlen *= treeinfo->brlen_scalers[p];

        int ret = pll_update_prob_matrices (treeinfo->partitions[p],
                                  treeinfo->param_indices[p],
                                  &m,
                                  &p_brlen,
                                  1);

        if (!ret)
          return PLL_FAILURE;

        treeinfo->pmatrix_valid[p][m] = 1;
        updated++;
      }
    }
  }

  return PLL_SUCCESS;
}

PLL_EXPORT void pllmod_treeinfo_invalidate_all(pllmod_treeinfo_t * treeinfo)
{
  unsigned int i, m;
  unsigned int clv_count = treeinfo->tip_count + (treeinfo->tip_count - 2) * 3;
  unsigned int pmatrix_count = treeinfo->tree->edge_count;

  for (i = 0; i < treeinfo->init_partition_count; ++i)
  {
    unsigned int p = treeinfo->init_partition_idx[i];

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
                                             pll_unode_t ** travbuffer,
                                             unsigned int travbuffer_size)
{
  for (unsigned int i = 0; i < treeinfo->init_partition_count; ++i)
  {
    unsigned int p = treeinfo->init_partition_idx[i];

    /* only selected partitioned will be affected */
    if (treeinfo_partition_active(treeinfo, p))
    {
      for (unsigned int j = 0; j < travbuffer_size; ++j)
      {
        const pll_unode_t * node = travbuffer[j];
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
                                                   const pll_unode_t * edge)
{
  for (unsigned int i = 0; i < treeinfo->init_partition_count; ++i)
  {
    unsigned int p = treeinfo->init_partition_idx[i];
    if (treeinfo->pmatrix_valid[p] && treeinfo_partition_active(treeinfo, p))
      treeinfo->pmatrix_valid[p][edge->pmatrix_index] = 0;
  }
}

PLL_EXPORT void pllmod_treeinfo_invalidate_clv(pllmod_treeinfo_t * treeinfo,
                                               const pll_unode_t * edge)
{
  for (unsigned int i = 0; i < treeinfo->init_partition_count; ++i)
  {
    unsigned int p = treeinfo->init_partition_idx[i];
    if (treeinfo->clv_valid[p] && treeinfo_partition_active(treeinfo, p))
      treeinfo->clv_valid[p][edge->node_index] = 0;
  }
}

static double treeinfo_compute_loglh(pllmod_treeinfo_t * treeinfo,
                                     int incremental,
                                     int update_pmatrices,
                                     double ** persite_lnl)
{
  /* tree root must be an inner node! */
  assert(!pllmod_utree_is_tip(treeinfo->root));

  unsigned int traversal_size = 0;
  unsigned int ops_count;
  unsigned int i, p;

  const double LOGLH_NONE = (double) NAN;
  double total_loglh = 0.0;
  const int old_active_partition = treeinfo->active_partition;

  /* NOTE: in unlinked brlen mode, up-to-date brlens for partition p
   * have to be prefetched to treeinfo->branch_lengths[p] !!! */
  int collect_brlen =
      (treeinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED ? 0 : 1);

  pllmod_treeinfo_set_active_partition(treeinfo, PLLMOD_TREEINFO_PARTITION_ALL);

  /* we need full traversal in 2 cases: 1) update p-matrices, 2) update all CLVs */
  if (!incremental || (update_pmatrices && collect_brlen))
  {
    /* perform a FULL postorder traversal of the unrooted tree */
    if (!pll_utree_traverse(treeinfo->root,
                            PLL_TREE_TRAVERSE_POSTORDER,
                            cb_full_traversal,
                            treeinfo->travbuffer,
                            &traversal_size))
      return LOGLH_NONE;
  }

  /* update p-matrices if asked for */
  if (update_pmatrices)
  {
    if (collect_brlen)
    {
      assert(traversal_size == treeinfo->tip_count * 2 - 2);
      for (i = 0; i < traversal_size; ++i)
       {
         pll_unode_t * node = treeinfo->travbuffer[i];
         treeinfo->branch_lengths[0][node->pmatrix_index] = node->length;
       }
    }

    pllmod_treeinfo_update_prob_matrices(treeinfo, !incremental);
  }

  if (incremental)
  {
    /* compute partial traversal and update only invalid CLVs */
    if (!pll_utree_traverse(treeinfo->root,
                            PLL_TREE_TRAVERSE_POSTORDER,
                            cb_partial_traversal,
                            treeinfo->travbuffer,
                            &traversal_size))
      return LOGLH_NONE;
  }

  /* create operations based on partial traversal obtained above */
  pll_utree_create_operations(treeinfo->travbuffer,
                              traversal_size,
                              NULL,
                              NULL,
                              treeinfo->operations,
                              NULL,
                              &ops_count);

  treeinfo->counter += ops_count;

//  printf("Traversal size (%s): %u\n", incremental ? "part" : "full", ops_count);

  /* iterate over all partitions (we assume that traversal is the same) */
  for (p = 0; p < treeinfo->partition_count; ++p)
  {
    if (!treeinfo->partitions[p])
    {
      /* this partition will be computed by another thread(s) */
      treeinfo->partition_loglh[p] = 0.0;
      continue;
    }

    /* all subsequent operation will affect current partition only */
    pllmod_treeinfo_set_active_partition(treeinfo, (int)p);

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
    treeinfo->partition_loglh[p] = pll_compute_edge_loglikelihood(
                                            treeinfo->partitions[p],
                                            treeinfo->root->clv_index,
                                            treeinfo->root->scaler_index,
                                            treeinfo->root->back->clv_index,
                                            treeinfo->root->back->scaler_index,
                                            treeinfo->root->pmatrix_index,
                                            treeinfo->param_indices[p],
                                            persite_lnl ? persite_lnl[p] : NULL);
  }

  /* sum up likelihood from all threads */
  if (treeinfo->parallel_reduce_cb)
  {
    treeinfo->parallel_reduce_cb(treeinfo->parallel_context,
                                 treeinfo->partition_loglh,
                                 p,
                                 PLLMOD_COMMON_REDUCE_SUM);
  }

  /* accumulate loglh by summing up over all the partitions */
  for (p = 0; p < treeinfo->partition_count; ++p)
    total_loglh += treeinfo->partition_loglh[p];

  /* restore original active partition */
  pllmod_treeinfo_set_active_partition(treeinfo, old_active_partition);

  assert(total_loglh < 0.);

  return total_loglh;
}

PLL_EXPORT double pllmod_treeinfo_compute_loglh(pllmod_treeinfo_t * treeinfo,
                                                int incremental)
{
  return treeinfo_compute_loglh(treeinfo, incremental, 1, NULL);
}

PLL_EXPORT double pllmod_treeinfo_compute_loglh_flex(pllmod_treeinfo_t * treeinfo,
                                                     int incremental,
                                                     int update_pmatrices)
{
  return treeinfo_compute_loglh(treeinfo, incremental, update_pmatrices, NULL);
}

PLL_EXPORT double pllmod_treeinfo_compute_loglh_persite(pllmod_treeinfo_t * treeinfo,
                                                        int incremental,
                                                        double ** persite_lnl)
{
  return treeinfo_compute_loglh(treeinfo, incremental, 1, persite_lnl);
}

PLL_EXPORT
int pllmod_treeinfo_scale_branches_all(pllmod_treeinfo_t * treeinfo, double scaler)
{
  unsigned int i, p;

  for (i = 0; i < treeinfo->subnode_count; ++i)
    treeinfo->subnodes[i]->length *= scaler;

  if (treeinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED)
  {
    for (p = 0; p < treeinfo->partition_count; ++p)
    {
      for (i = 0; i < treeinfo->tree->edge_count; ++i)
        treeinfo->branch_lengths[p][i] *= scaler;
    }

  }
  else
  {
    for (i = 0; i < treeinfo->tree->edge_count; ++i)
      treeinfo->branch_lengths[0][i] *= scaler;
  }

  return PLL_SUCCESS;
}

PLL_EXPORT
int pllmod_treeinfo_scale_branches_partition(pllmod_treeinfo_t * treeinfo,
                                             unsigned int partition_idx,
                                             double scaler)
{
  unsigned int i;

  if (treeinfo->brlen_linkage != PLLMOD_COMMON_BRLEN_UNLINKED)
  {
    pllmod_set_error(PLLMOD_TREE_ERROR_INVALID_TREE,
                     "Per-partition branch length scaling is supported "
                     "in unlinked branch length mode only.\n");
    return PLL_FAILURE;
  }

  if (partition_idx >= treeinfo->partition_count)
  {
    pllmod_set_error(PLL_ERROR_PARAM_INVALID,
                     "Partition ID out of bounds\n");
    return PLL_FAILURE;
  }

  for (i = 0; i < treeinfo->tree->edge_count; ++i)
    treeinfo->branch_lengths[partition_idx][i] *= scaler;

  return PLL_SUCCESS;
}

PLL_EXPORT
int pllmod_treeinfo_normalize_brlen_scalers(pllmod_treeinfo_t * treeinfo)
{
  double sum_scalers = 0.;
  double sum_sites = 0.;
  unsigned int p;

  for (p = 0; p < treeinfo->partition_count; ++p)
  {
    if (treeinfo->partitions[p])
    {
      const double pat_sites = treeinfo->partitions[p]->pattern_weight_sum;
      sum_sites += pat_sites;
      sum_scalers += treeinfo->brlen_scalers[p] * pat_sites;
    }
  }

  /* sum up scalers and sites from all threads */
  if (treeinfo->parallel_reduce_cb)
  {
    // TODO: although this reduce is unlikely to become a bottleneck, it can
    // be avoided by storing _total_ pattern_weight_sum in treeinfo
    treeinfo->parallel_reduce_cb(treeinfo->parallel_context,
                                 &sum_scalers,
                                 1,
                                 PLLMOD_COMMON_REDUCE_SUM);

    treeinfo->parallel_reduce_cb(treeinfo->parallel_context,
                                 &sum_sites,
                                 1,
                                 PLLMOD_COMMON_REDUCE_SUM);
  }

  const double mean_rate = sum_scalers / sum_sites;
  pllmod_treeinfo_scale_branches_all(treeinfo, mean_rate);
  for (p = 0; p < treeinfo->partition_count; ++p)
  {
    if (treeinfo->partitions[p])
      treeinfo->brlen_scalers[p] /= mean_rate;
  }

  return PLL_SUCCESS;
}

static int treeinfo_check_tree(pllmod_treeinfo_t * treeinfo,
                               pll_utree_t * tree)
{
  if (!treeinfo || !tree)
  {
    pllmod_set_error(PLL_ERROR_PARAM_INVALID,
                     "Parameter is NULL\n");
    return PLL_FAILURE;
  }

  if (treeinfo->tip_count != tree->tip_count)
  {
    pllmod_set_error(PLLMOD_TREE_ERROR_INVALID_TREE_SIZE,
                     "Invalid tree size. Got %d instead of %d\n",
                     tree->tip_count, treeinfo->tip_count);
    return PLL_FAILURE;
  }

  if (!tree->binary)
  {
    pllmod_set_error(PLLMOD_TREE_ERROR_INVALID_TREE,
                     "Binary tree expected\n");
    return PLL_FAILURE;
  }

  return PLL_SUCCESS;
}

static int treeinfo_init_tree(pllmod_treeinfo_t * treeinfo)
{
  pll_utree_t * tree = treeinfo->tree;
  unsigned int node_count = tree->tip_count + tree->inner_count;
  unsigned int edge_count = tree->edge_count;

  /* save virtual root */
  treeinfo->root = tree->vroot;

  /* 1. save back pointer to treeinfo stuct in in eacn node's _data_ field
   * 2. collect branch length from the tree and store in a separate array
   *    indexed by pmatrix_index */
  for (unsigned int i = 0; i < node_count; ++i)
  {
    pll_unode_t * snode = tree->nodes[i];
    do
    {
      unsigned int pmat_idx = snode->pmatrix_index;
      if (pmat_idx > edge_count)
      {
        pllmod_set_error(PLLMOD_TREE_ERROR_INVALID_TREE,
                         "p-matrix index out of bounds (%u). "
                         "treeinfo structure require that each branch "
                         "branch is assigned a unique p-matrix index "
                         "between 0 and branch_count-1\n",
                         pmat_idx);
        return PLL_FAILURE;
      }
      treeinfo->subnodes[snode->node_index] = snode;
      treeinfo->branch_lengths[0][pmat_idx] = snode->length;
      snode->data = treeinfo;
      snode = snode->next;
    }
    while (snode && snode != tree->nodes[i]);
  }

  /* in unlinked branch length mode, copy brlen to other partitions */
  if (treeinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED)
  {
    for (unsigned int i = 1; i < treeinfo->partition_count; ++i)
    {
      // TODO: only save brlens for initialized partitions
//      if (treeinfo->partitions[i])
      {
        memcpy(treeinfo->branch_lengths[i], treeinfo->branch_lengths[0],
               edge_count * sizeof(double));
      }
    }
  }

  return PLL_SUCCESS;
}

PLL_EXPORT int pllmod_treeinfo_set_tree(pllmod_treeinfo_t * treeinfo,
                                        pll_utree_t * tree)
{
  if (!treeinfo_check_tree(treeinfo, tree))
    return PLL_FAILURE;

  unsigned int node_count = tree->tip_count + tree->inner_count;

  /* make a shallow copy of pll_utree_t structure, that is,
   * pll_unodes will not be cloned! */
  pll_unode_t ** nodes = treeinfo->tree->nodes;
  memcpy(treeinfo->tree, tree, sizeof(pll_utree_t));
  treeinfo->tree->nodes = nodes;
  memcpy(treeinfo->tree->nodes, tree->nodes, node_count*sizeof(pll_unode_t *));

  return treeinfo_init_tree(treeinfo);
}

PLL_EXPORT int pllmod_treeinfo_set_constraint_clvmap(pllmod_treeinfo_t * treeinfo,
                                                     const int * clv_index_map)
{
  const unsigned int tip_count = treeinfo->tree->tip_count;
  const unsigned int inner_count = treeinfo->tree->inner_count;
  const size_t cons_size = (tip_count + inner_count) * sizeof(unsigned int);

  if(!treeinfo->constraint)
  {
    treeinfo->constraint = malloc(cons_size);
    if (!treeinfo->constraint)
    {
      pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                       "Can't allocate memory for topological constraint\n");
      return PLL_FAILURE;
    }
  }

  assert(treeinfo->constraint);

  for (unsigned int i = 0; i < tip_count + inner_count; ++i)
  {
    const pll_unode_t * node = treeinfo->tree->nodes[i];
    const unsigned int cons_group_id = (unsigned int) clv_index_map[node->clv_index]+1;

    treeinfo->constraint[node->clv_index] = cons_group_id;
  }

  return PLL_SUCCESS;
}

PLL_EXPORT int pllmod_treeinfo_set_constraint_tree(pllmod_treeinfo_t * treeinfo,
                                                   const pll_utree_t * cons_tree)
{
  unsigned int node_count = cons_tree->tip_count * 2 - 2;
  int * clv_index_map = NULL;
  int retval;

  if (treeinfo->tip_count < cons_tree->tip_count)
  {
    pllmod_set_error(PLLMOD_TREE_ERROR_INVALID_TREE_SIZE,
                     "Invalid tree size. Got %d instead of %d\n",
                     cons_tree->tip_count, treeinfo->tip_count);
    return PLL_FAILURE;
  }

  clv_index_map = (int *) calloc(node_count, sizeof(int));
  if (!clv_index_map)
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Can't allocate memory for clv_index_map\n");
    return PLL_FAILURE;
  }

  pll_utree_t * bin_cons_tree = pllmod_utree_resolve_multi(cons_tree,
                                                           0, clv_index_map);

  if (!bin_cons_tree)
  {
    assert(pll_errno);
    free(clv_index_map);
    return PLL_FAILURE;
  }

  // HACK; copy clv group ids from inner nodes to tips
  for (unsigned int i = 0; i < bin_cons_tree->tip_count; ++i)
  {
    pll_unode_t * node = bin_cons_tree->nodes[i];
    assert(!node->next);
    clv_index_map[node->clv_index] = clv_index_map[node->back->clv_index];
  }

  pll_utree_destroy(bin_cons_tree, NULL);

  if (treeinfo->tip_count > cons_tree->tip_count)
  {
    // non-comprehensive constraint
    unsigned int free_count = treeinfo->tip_count - cons_tree->tip_count;
    unsigned int ext_node_count = treeinfo->tip_count * 2 - 2;
    int * ext_clv_index_map = (int *) calloc(ext_node_count, sizeof(int));

    for (unsigned int i = 0; i < ext_node_count; ++i)
    {
      unsigned int clv_id = treeinfo->tree->nodes[i]->clv_index;
      if (clv_id < cons_tree->tip_count)
      {
        // copy map for tips in the constraint tree and adjust clv_ids
        ext_clv_index_map[clv_id] = clv_index_map[clv_id] + free_count;
      }
      else if (clv_id < treeinfo->tip_count)
      {
        // mark new tips as freely movable
        ext_clv_index_map[clv_id] = -1;
      }
      else if (clv_id < node_count + free_count)
      {
        // update clv_ids for old inner nodes by adding free tip count
        unsigned int old_clv_id = clv_id - free_count;
        ext_clv_index_map[clv_id] = clv_index_map[old_clv_id] + free_count;
      }
      else if (clv_id < ext_node_count)
      {
        // mark new inner nodes as freely movable
        ext_clv_index_map[clv_id] = -1;
      }
      else
        assert(0);
    }

    free(clv_index_map);
    clv_index_map = ext_clv_index_map;
  }

  retval = pllmod_treeinfo_set_constraint_clvmap(treeinfo, clv_index_map);

  free(clv_index_map);

  return retval;
}

static unsigned int find_cons_id(pll_unode_t * node,
                                 const unsigned int * constraint,
                                 unsigned int s)
{
  unsigned int cons_group_id = constraint[node->clv_index];
  if (!node->next || cons_group_id > 0)
    return cons_group_id;
  else
  {
    unsigned int left_id = find_cons_id(node->next->back, constraint, s);
    unsigned int right_id = find_cons_id(node->next->next->back, constraint, s);
    if (left_id == right_id)
      return left_id;
    else if (!s)
      return PLL_MAX(left_id, right_id);
    else
      return (left_id == 0 || left_id == s) ? right_id : left_id;
  }
}

PLL_EXPORT int pllmod_treeinfo_check_constraint(pllmod_treeinfo_t * treeinfo,
                                                pll_unode_t * subtree,
                                                pll_unode_t * regraft_edge)
{
  if (treeinfo->constraint)
  {
    int res;
    unsigned int s = treeinfo->constraint[subtree->clv_index];
    s  = s ? s : find_cons_id(subtree->back, treeinfo->constraint, 0);

    if (s)
    {
      unsigned int r1 = find_cons_id(regraft_edge, treeinfo->constraint, s);
      unsigned int r2 = find_cons_id(regraft_edge->back, treeinfo->constraint, s);

      res = (s == r1 || s == r2) ? PLL_SUCCESS : PLL_FAILURE;
    }
    else
      res = PLL_SUCCESS;

    return res;
  }
  else
    return PLL_SUCCESS;
}


static pllmod_ancestral_t * pllmod_treeinfo_create_ancestral(const pllmod_treeinfo_t * treeinfo)
{
  unsigned int i;

  pllmod_ancestral_t * ancestral = (pllmod_ancestral_t *) calloc(1, sizeof(pllmod_ancestral_t));

  if (!ancestral)
    return NULL;

  ancestral->node_count = treeinfo->tree->inner_count;
  ancestral->partition_count = treeinfo->init_partition_count;

  ancestral->nodes = (pll_unode_t **) calloc(ancestral->node_count,
                                             sizeof(pll_unode_t *));
  ancestral->partition_indices = (unsigned int *) calloc(ancestral->partition_count,
                                                         sizeof(unsigned int));

  ancestral->probs = (double **) calloc(ancestral->node_count, sizeof(double *));

  size_t vector_size = 0;
  for (i = 0; i < ancestral->partition_count; ++i)
  {
    const pll_partition_t * partition = treeinfo->init_partitions[i];
    vector_size += partition->sites * partition->states;
  }


  ancestral->tree = pll_utree_clone(treeinfo->tree);

  for (i = 0; i < ancestral->node_count; ++i)
    ancestral->probs[i] = (double *) calloc(vector_size, sizeof(double));

  return ancestral;
}

PLL_EXPORT void pllmod_treeinfo_destroy_ancestral(pllmod_ancestral_t * ancestral)
{
  unsigned int i;

  if (!ancestral)
    return;

  free(ancestral->nodes);
  free(ancestral->partition_indices);

  for (i = 0; i < ancestral->node_count; ++i)
    free(ancestral->probs[i]);

  pll_utree_destroy(ancestral->tree, NULL);
  free(ancestral->probs);
}

PLL_EXPORT
pllmod_ancestral_t * pllmod_treeinfo_compute_ancestral(pllmod_treeinfo_t * treeinfo)
{
  unsigned int i, p;
  unsigned int traversal_size;
  unsigned int node_count = treeinfo->tree->tip_count + treeinfo->tree->inner_count;

  pll_unode_t * old_root = treeinfo->root;

  pll_unode_t ** travbuffer = (pll_unode_t **) calloc(node_count, sizeof(pll_unode_t *));

  pllmod_ancestral_t * ancestral = pllmod_treeinfo_create_ancestral(treeinfo);

  if (!travbuffer || !ancestral)
  {
    if (travbuffer)
      free(travbuffer);
    if (ancestral)
      pllmod_treeinfo_destroy_ancestral(ancestral);
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Can't allocate memory for ancestral probabilities\n");
    return NULL;
  }

  /* perform a FULL postorder traversal of the unrooted tree */
  if (!pll_utree_traverse(ancestral->tree->vroot,
                          PLL_TREE_TRAVERSE_POSTORDER,
                          cb_full_traversal,
                          travbuffer,
                          &traversal_size))
  {
    free(travbuffer);
    pllmod_treeinfo_destroy_ancestral(ancestral);
    return NULL;
  }

  assert(traversal_size == node_count);

#ifdef DEBUG
  printf("\nTraversal: ");
  for (i = 0; i < traversal_size; ++i)
    printf("%u ", travbuffer[i]->clv_index);
  printf("\n\n");
#endif

  /* collect internal nodes from the traversal */
  unsigned int n = 0;
  for (i = 0; i < traversal_size; ++i)
  {
    pll_unode_t * node = travbuffer[i];

    if (pllmod_utree_is_tip(node))
      continue;

    /* generate inner node label if it is empty */
    if (!node->label)
    {
      char buf[100];
      snprintf(buf, 100, "Node%u", n+1);
      node->label = node->next->label = node->next->next->label = strdup(buf);
    }

    ancestral->nodes[n] = node;
    ++n;
  }

  assert(n == ancestral->node_count);

  free(travbuffer);

  /* now compute ancestral state probs for every internal node */
  for (i = 0; i < ancestral->node_count; ++i)
  {
    pll_unode_t * node = ancestral->nodes[i];
    pll_unode_t * treeinfo_node = treeinfo->subnodes[node->node_index];
    double * ancp = ancestral->probs[i];

    treeinfo->root = treeinfo_node;
    pllmod_treeinfo_compute_loglh(treeinfo, 1);

    for (p = 0; p < treeinfo->init_partition_count; ++p)
    {
      pll_partition_t * partition = treeinfo->init_partitions[p];
      size_t part_span = partition->sites * partition->states;
      unsigned int pidx = treeinfo->init_partition_idx[p];

      ancestral->partition_indices[p] = pidx;

      if (!pll_compute_node_ancestral(partition,
                                      node->clv_index,
                                      node->scaler_index,
                                      node->back->clv_index,
                                      node->back->scaler_index,
                                      node->pmatrix_index,
                                      treeinfo->param_indices[pidx],
                                      ancp))
      {
        pllmod_treeinfo_destroy_ancestral(ancestral);
        return NULL;
      }

      ancp += part_span;
    }
  }

  treeinfo->root = old_root;

  return ancestral;
}
