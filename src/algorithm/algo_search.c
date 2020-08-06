/*
Copyright (C) 2016 Alexey Kozlov, Diego Darriba, Tomas Flouri and Alexandros Stamatakis.

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

 /**
  * @file algo_search.c
  *
  * @brief Core topological search algorithms
  *
  * @author Alexey Kozlov
  * @author Diego Darriba
  */

#include "pllmod_algorithm.h"
#include "../pllmod_common.h"

/* if not defined, branch length optimization will use
* the same starting set of branch lengths for every topology */
#define PLLMOD_SEARCH_GREEDY_BLO

/* if defined, branch length array will be always dynamically allocated,
* i.e. also in linked and scaled brlen modes; otherwise only in unlinked mode */
//#define PLLMOD_SEARCH_BRLEN_DYNALLOC

#define BRLEN_BUF_COUNT 12

typedef struct spr_params
{
  pll_bool_t thorough;
  unsigned int radius_min;
  unsigned int radius_max;
  unsigned int ntopol_keep;
  double bl_min;
  double bl_max;
  int smoothings;
  int brlen_opt_method;
  double * brlen_buf[BRLEN_BUF_COUNT];
} pllmod_search_params_t;

typedef struct rollback_list
{
  pll_tree_rollback_t * list;
  size_t current;
  unsigned int round;
  size_t size;
} pllmod_rollback_list_t;

typedef struct node_entry {
  pll_unode_t * p_node;
  pll_unode_t * r_node;
  double bb1, bb2, bb3;
  double *b1, *b2, *b3;
  double lh;
  unsigned int rollback_num;
} node_entry_t;

typedef struct bestnode_list
{
  node_entry_t * list;
  size_t current;
  size_t size;
  unsigned int brlen_set_count;
  double ** brlen_buffers;
} pllmod_bestnode_list_t;

static void algo_query_allnodes_recursive(pll_unode_t * node,
                                          pll_unode_t ** buffer,
                                          unsigned int * index)
{
  if (node->next)
  {
    algo_query_allnodes_recursive(node->next->back, buffer, index);
    algo_query_allnodes_recursive(node->next->next->back, buffer, index);

    buffer[(*index)++] = node->next->next;
    buffer[(*index)++] = node->next;
    buffer[(*index)++] = node;
  }
}

static unsigned int algo_query_allnodes(pll_unode_t * root, pll_unode_t ** buffer)
{
  assert(root && buffer);

  unsigned int index = 0;

  algo_query_allnodes_recursive(root->back, buffer, &index);
  algo_query_allnodes_recursive(root, buffer, &index);

  return index;
}

static pllmod_rollback_list_t * algo_rollback_list_create(size_t slots)
{
  pllmod_rollback_list_t * rollback_list =
    (pllmod_rollback_list_t *) calloc(1, sizeof(pllmod_rollback_list_t));
  if (!rollback_list)
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for rollback list\n");
    return NULL;
  }
  rollback_list->current = 0;
  rollback_list->round = 0;
  rollback_list->size = slots;
  if (slots > 0)
  {
    rollback_list->list =
      (pll_tree_rollback_t *) calloc(slots, sizeof(pll_tree_rollback_t));
    if (!rollback_list->list)
    {
      pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                       "Cannot allocate memory for rollback list items\n");
      free(rollback_list);
      return NULL;
    }
  }
  return rollback_list;
}

static void algo_rollback_list_destroy(pllmod_rollback_list_t * rollback_list)
{
  if (rollback_list)
  {
    if(rollback_list->list)
      free(rollback_list->list);
    free(rollback_list);
  }
}

static pll_tree_rollback_t * algo_rollback_list_prev(
                               pllmod_rollback_list_t * rollback_list)
{
  if (rollback_list->current > 0)
    rollback_list->current--;
  else if (rollback_list->round > 0 && rollback_list->size > 0)
  {
    rollback_list->round--;
    rollback_list->current = rollback_list->size ? rollback_list->size - 1 : 0;
  }
  else
    return NULL;

  return rollback_list->list + rollback_list->current;
}

static pll_tree_rollback_t * algo_rollback_list_next(
                               pllmod_rollback_list_t * rollback_list)
{
  if (rollback_list->current + 1 < rollback_list->size)
    rollback_list->current++;
  else
  {
    rollback_list->round++;
    rollback_list->current = 0;
  }

  return rollback_list->list + rollback_list->current;
}

static int algo_rollback_list_abspos(pllmod_rollback_list_t * rollback_list)
{
  return rollback_list->size * rollback_list->round + rollback_list->current;
}

/*                  *
 *  best_node_list  *
 *                  */

static void algo_bestnode_list_destroy(pllmod_bestnode_list_t * bestnode_list)
{
 if (bestnode_list)
 {
   if (bestnode_list->brlen_buffers)
   {
     for (size_t i = 0; i < bestnode_list->size; ++i)
       free(bestnode_list->brlen_buffers[i]);
     free(bestnode_list->brlen_buffers);
   }
   if (bestnode_list->list)
     free(bestnode_list->list);
   free(bestnode_list);
 }
}

static pllmod_bestnode_list_t * algo_bestnode_list_create(size_t slots,
                                                          unsigned int brlen_set_count)
{
  pllmod_bestnode_list_t * bestnode_list =
    (pllmod_bestnode_list_t *) calloc(1, sizeof(pllmod_bestnode_list_t));
  if (!bestnode_list)
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for best node list\n");
    return NULL;
  }
  bestnode_list->current = 0;
  bestnode_list->size = slots;
  bestnode_list->brlen_set_count = brlen_set_count;
  if (slots > 0)
  {
    bestnode_list->list = (node_entry_t *) calloc(slots, sizeof(node_entry_t));
    if (!bestnode_list->list)
    {
      pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                       "Cannot allocate memory for best node list items\n");
      free(bestnode_list);
      return NULL;
    }

#ifndef PLLMOD_SEARCH_BRLEN_DYNALLOC
    if (brlen_set_count == 1)
    {
      bestnode_list->brlen_buffers = NULL;
      for (size_t i = 0; i < slots; ++i)
      {
        bestnode_list->list[i].b1 = &bestnode_list->list[i].bb1;
        bestnode_list->list[i].b2 = &bestnode_list->list[i].bb2;
        bestnode_list->list[i].b3 = &bestnode_list->list[i].bb3;
      }
    }
    else
#endif
    {
      /* allocate branch length buffers to store unlinked branches */
      bestnode_list->brlen_buffers = (double **) calloc(slots, sizeof(double *));
      for (size_t i = 0; i < slots; ++i)
      {
        bestnode_list->brlen_buffers[i] = (double *) calloc(3*brlen_set_count,
                                                            sizeof(double));
      }
      if (!bestnode_list->brlen_buffers || !bestnode_list->brlen_buffers[slots-1])
      {
        pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                         "Cannot allocate memory for best node list items\n");
        algo_bestnode_list_destroy(bestnode_list);
        return NULL;
      }

      /* set branch lengths pointers to the corresponding slots */
      for (size_t i = 0; i < slots; ++i)
      {
        bestnode_list->list[i].b1 = bestnode_list->brlen_buffers[i];
        bestnode_list->list[i].b2 =
            bestnode_list->brlen_buffers[i] + brlen_set_count;
        bestnode_list->list[i].b3 =
            bestnode_list->brlen_buffers[i] + 2*brlen_set_count;
      }
    }
  }
  return bestnode_list;
}

static void algo_bestnode_list_copy_entry(pllmod_bestnode_list_t * best_node_list,
                                          size_t idx,
                                          const node_entry_t * src)
{
  node_entry_t * dst = &best_node_list->list[idx];
  dst->p_node = src->p_node;
  dst->r_node = src->r_node;
  dst->lh = src->lh;
  dst->rollback_num = src->rollback_num;

  if (best_node_list->brlen_buffers)
  {
    memcpy(dst->b1, src->b1, best_node_list->brlen_set_count * sizeof(double));
    memcpy(dst->b2, src->b2, best_node_list->brlen_set_count * sizeof(double));
    memcpy(dst->b3, src->b3, best_node_list->brlen_set_count * sizeof(double));
  }
  else
  {
    dst->bb1 = *(src->b1);
    dst->bb2 = *(src->b2);
    dst->bb3 = *(src->b3);
  }
}

static void algo_bestnode_list_save(pllmod_bestnode_list_t * best_node_list,
                                    const node_entry_t * entry)
{
  node_entry_t * list = best_node_list->list;
  const size_t list_size = best_node_list->size;
  size_t idx = 0, j;

  while (idx < list_size && (list[idx].p_node) && (entry->lh < list[idx].lh))
  {
    ++idx;
  }

  /* do not insert: candidate tree LH too low, or list has size of 0 */
  if (idx >= list_size)
    return;

  j = list_size - 1;
  while (j > idx)
  {
    algo_bestnode_list_copy_entry(best_node_list, j, &list[j-1]);
    --j;
  }

  assert(idx >= 0 && idx < list_size);

  algo_bestnode_list_copy_entry(best_node_list, idx, entry);
}

static int algo_bestnode_list_next_index(pllmod_bestnode_list_t * best_node_list,
                                         unsigned int rollback_num,
                                         int curr_index)
{
  assert(curr_index >= -1);

  node_entry_t * list = best_node_list->list;
  const unsigned int list_size = best_node_list->size;

  do
  {
    curr_index++;
    if (curr_index >= (int) list_size || !list[curr_index].p_node)
      return -1;
  }
  while (list[curr_index].rollback_num != rollback_num);

  return curr_index;
}

#ifdef DEBUG
static void algo_bestnode_list_print(pllmod_bestnode_list_t * best_node_list)
{
  printf("\nBESTNODE_LIST (%lu / %lu): \n", best_node_list->current, best_node_list->size);
  for (size_t i = 0; i < best_node_list->size; ++i)
  {
    node_entry_t * entry = &best_node_list->list[i];
    printf("%.9lf  %.9lf   %.9lf\n", entry->b1[0], entry->b2[0], entry->b3[0]);
  }
  printf("\n");
}
#endif


static double algo_optimize_bl_iterative(pll_unode_t * node,
                                         pllmod_treeinfo_t * treeinfo,
                                         const pllmod_search_params_t * params,
                                         int radius,
                                         double lh_epsilon,
                                         double smooth_factor)
{
  int smoothings = (int) round(smooth_factor * params->smoothings);

  double new_loglh = pllmod_opt_optimize_branch_lengths_local_multi(
                                                  treeinfo->partitions,
                                                  treeinfo->partition_count,
                                                  node,
                                                  treeinfo->param_indices,
                                                  treeinfo->deriv_precomp,
                                                  treeinfo->branch_lengths,
                                                  treeinfo->brlen_scalers,
                                                  params->bl_min,
                                                  params->bl_max,
                                                  lh_epsilon,
                                                  smoothings,
                                                  radius,
                                                  1,       /* keep_update */
                                                  params->brlen_opt_method,
                                                  treeinfo->brlen_linkage,
                                                  treeinfo->parallel_context,
                                                  treeinfo->parallel_reduce_cb);

  if (new_loglh)
    return -1 * new_loglh;
  else
  {
    assert(pll_errno);
    return 0;
  }
}

static double algo_optimize_bl_triplet(pll_unode_t * node,
                                       pllmod_treeinfo_t * treeinfo,
                                       const pllmod_search_params_t * params,
                                       double smooth_factor)
{
  return algo_optimize_bl_iterative(node, treeinfo, params,
                                    1, 0.1, smooth_factor);
}

static double algo_optimize_bl_all(pllmod_treeinfo_t * treeinfo,
                                   const pllmod_search_params_t * params,
                                   double lh_epsilon,
                                   double smooth_factor)
{
  pllmod_treeinfo_compute_loglh(treeinfo, 0);

  return algo_optimize_bl_iterative(treeinfo->root, treeinfo, params,
                                    PLLMOD_OPT_BRLEN_OPTIMIZE_ALL, lh_epsilon,
                                    smooth_factor);
}

static void algo_unode_fix_length(pllmod_treeinfo_t * treeinfo,
                                  pll_unode_t * node,
                                  double bl_min, double bl_max)
{
  unsigned int p;
  unsigned int pmatrix_index = node->pmatrix_index;

  if (treeinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED)
  {
    for (p = 0; p < treeinfo->partition_count; ++p)
    {
      if (treeinfo->partitions[p])
      {
        double p_brlen = treeinfo->branch_lengths[p][pmatrix_index];
        if (p_brlen < bl_min)
        {
          pllmod_treeinfo_set_branch_length_partition(treeinfo, node, p, bl_min);
          treeinfo->pmatrix_valid[p][pmatrix_index] = 0;
        }
        else if (p_brlen > bl_max)
        {
          pllmod_treeinfo_set_branch_length_partition(treeinfo, node, p, bl_max);
          treeinfo->pmatrix_valid[p][pmatrix_index] = 0;
        }
      }
    }
  }
  else
  {
    if (node->length < bl_min)
    {
      pllmod_treeinfo_set_branch_length(treeinfo, node, bl_min);
      pllmod_treeinfo_invalidate_pmatrix(treeinfo, node);
    }
    else if (node->length > bl_max)
    {
      pllmod_treeinfo_set_branch_length(treeinfo, node, bl_max);
      pllmod_treeinfo_invalidate_pmatrix(treeinfo, node);
    }
  }
}

int algo_update_pmatrix(pllmod_treeinfo_t * treeinfo,
                        pll_unode_t * edge)
{
  unsigned int p;
  unsigned int updated = 0;
  unsigned int pmatrix_index = edge->pmatrix_index;

  for (p = 0; p < treeinfo->partition_count; ++p)
  {
    /* only selected partitioned will be affected */
    if (treeinfo->partitions[p])
    {
      if (treeinfo->pmatrix_valid[p][pmatrix_index])
        continue;

      double p_brlen = treeinfo->branch_lengths[p][pmatrix_index];
      if (treeinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_SCALED)
        p_brlen *= treeinfo->brlen_scalers[p];

      int ret = pll_update_prob_matrices (treeinfo->partitions[p],
                                          treeinfo->param_indices[p],
                                          &pmatrix_index,
                                          &p_brlen,
                                          1);

      if (!ret)
        return PLL_FAILURE;

      treeinfo->pmatrix_valid[p][pmatrix_index] = 1;
      updated++;
    }
  }

  return PLL_SUCCESS;
}

static pll_unode_t * algo_utree_prune(pllmod_treeinfo_t * treeinfo,
                                      const pllmod_search_params_t * params,
                                      pll_unode_t * edge)
{
  pll_unode_t * orig_prune_edge;
  if (treeinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED)
  {
    double * b1 = params->brlen_buf[10];
    double * b2 = params->brlen_buf[11];

    pllmod_treeinfo_get_branch_length_all(treeinfo, edge->next->back, b1);
    pllmod_treeinfo_get_branch_length_all(treeinfo, edge->next->next->back, b2);

    orig_prune_edge = pllmod_utree_prune(edge);

    for (unsigned int i = 0; i < treeinfo->init_partition_count; ++i)
      b1[i] = b1[i] + b2[i];

    pllmod_treeinfo_set_branch_length_all(treeinfo, orig_prune_edge, b1);
  }
  else
  {
    double b1 = edge->next->next->back->length + edge->next->back->length;

    orig_prune_edge = pllmod_utree_prune(edge);

    pllmod_treeinfo_set_branch_length_all(treeinfo, orig_prune_edge, &b1);
  }

  return orig_prune_edge;
}


static int algo_utree_regraft(pllmod_treeinfo_t * treeinfo,
                              const pllmod_search_params_t * params,
                              pll_unode_t * p_edge,
                              pll_unode_t * r_edge)
{
  int retval;
  if (treeinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED)
  {
    double * b1 = params->brlen_buf[10];

    pllmod_treeinfo_get_branch_length_all(treeinfo, r_edge, b1);

    for (unsigned int i = 0; i < treeinfo->init_partition_count; ++i)
      b1[i] = b1[i] / 2.;

    retval = pllmod_utree_regraft(p_edge, r_edge);

    pllmod_treeinfo_set_branch_length_all(treeinfo, p_edge->next->back, b1);
    pllmod_treeinfo_set_branch_length_all(treeinfo, p_edge->next->next->back, b1);
  }
  else
  {
    double b1 = r_edge->length;

    b1 = b1 / 2.;

    retval = pllmod_utree_regraft(p_edge, r_edge);

    pllmod_treeinfo_set_branch_length_all(treeinfo, p_edge->next->back, &b1);
    pllmod_treeinfo_set_branch_length_all(treeinfo, p_edge->next->next->back, &b1);
  }

  return retval;
}

PLL_EXPORT int algo_utree_spr(pllmod_treeinfo_t * treeinfo,
                              const pllmod_search_params_t * params,
                              pll_unode_t * p_edge,
                              pll_unode_t * r_edge,
                              pll_tree_rollback_t * rollback_info)
{
  int retval;

  if (pllmod_utree_is_tip(p_edge))
  {
    /* invalid move */
    pllmod_set_error(PLLMOD_TREE_ERROR_SPR_INVALID_NODE,
                     "Attempting to prune a leaf branch");
    return PLL_FAILURE;
  }

  // TODO do we need to store unlinked brlens here?
  /* save rollback information */
  if (rollback_info)
  {
    rollback_info->rearrange_type     = PLLMOD_TREE_REARRANGE_SPR;
    rollback_info->rooted             = 0;
    rollback_info->SPR.prune_edge     = (void *) p_edge;
    rollback_info->SPR.regraft_edge   = (void *) p_edge->next->back;
    rollback_info->SPR.prune_bl       = p_edge->length;
    rollback_info->SPR.prune_left_bl  = p_edge->next->length;
    rollback_info->SPR.prune_right_bl = p_edge->next->next->length;
    rollback_info->SPR.regraft_bl     = r_edge->length;
  }

  algo_utree_prune(treeinfo, params, p_edge);

  retval = algo_utree_regraft(treeinfo, params, p_edge, r_edge);

  return retval;
}

static int best_reinsert_edge(pllmod_treeinfo_t * treeinfo,
                              node_entry_t * entry,
                              cutoff_info_t * cutoff_info,
                              const pllmod_search_params_t * params)
{
  assert(treeinfo && entry && params);

  unsigned int i, j;
  pll_unode_t * orig_prune_edge;
  pll_unode_t ** regraft_nodes;
  pll_unode_t * r_edge;
  int regraft_edges;
  unsigned int r_dist;
  double *z1, *z2, *z3;
  double *b1, *b2, *b3;
  double * regraft_length;
  unsigned int redge_count = 0;
  unsigned int ncount;
  int retval;
  unsigned int * regraft_dist;
  int descent;
  double loglh;

  pll_unode_t * p_edge = entry->p_node;
  const size_t total_edge_count = treeinfo->tree->edge_count;

  entry->r_node = NULL;
  entry->lh = PLLMOD_OPT_LNL_UNLIKELY;

  /* init brlen buffer pointers */
  z1 = params->brlen_buf[0];
  z2 = params->brlen_buf[1];
  z3 = params->brlen_buf[2];

  b1 = params->brlen_buf[3];
  b2 = params->brlen_buf[4];
  b3 = params->brlen_buf[5];
  regraft_length = params->brlen_buf[6];

  /* save original branch lengths at the pruning point */
  pllmod_treeinfo_get_branch_length_all(treeinfo, p_edge, z1);
  pllmod_treeinfo_get_branch_length_all(treeinfo, p_edge->next, z2);
  pllmod_treeinfo_get_branch_length_all(treeinfo, p_edge->next->next, z3);

  pllmod_treeinfo_set_root(treeinfo, p_edge);

  /* recompute all CLVs and p-matrices before pruning */
  loglh = pllmod_treeinfo_compute_loglh(treeinfo, 0);

  /* PRUNE */
  orig_prune_edge = algo_utree_prune(treeinfo, params, p_edge);
  if (!orig_prune_edge)
  {
    /* check that errno was set correctly */
    assert(pll_errno & PLLMOD_TREE_ERROR_SPR_MASK);
    return PLL_FAILURE;
  }

  algo_unode_fix_length(treeinfo, orig_prune_edge, params->bl_min, params->bl_max);

  pllmod_treeinfo_set_root(treeinfo, orig_prune_edge);

  /* invalidate CLVs & p-matrix at the pruned edge */
  pllmod_treeinfo_invalidate_clv(treeinfo, orig_prune_edge);
  pllmod_treeinfo_invalidate_clv(treeinfo, orig_prune_edge->back);
  pllmod_treeinfo_invalidate_pmatrix(treeinfo, orig_prune_edge);

  /* recompute p-matrix for the original prune edge */
  algo_update_pmatrix(treeinfo, orig_prune_edge);

  /* get list of candidate regrafting nodes in the given distance range */
  regraft_nodes = (pll_unode_t **) calloc (total_edge_count,
                                           sizeof(pll_unode_t *));
  if (!regraft_nodes)
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for regraft nodes\n");
    return PLL_FAILURE;
  }

  retval = pllmod_utree_nodes_at_node_dist(treeinfo->root,
                                           &regraft_nodes[redge_count],
                                           &ncount,
                                           params->radius_min,
                                           params->radius_min);
  redge_count += ncount;

  if (!pllmod_utree_is_tip(treeinfo->root->back))
  {
    retval &= pllmod_utree_nodes_at_node_dist(treeinfo->root->back,
                                              &regraft_nodes[redge_count],
                                              &ncount,
                                              params->radius_min,
                                              params->radius_min);
    redge_count += ncount;
  }
  assert(retval == PLL_SUCCESS);

  /* initialize regraft distances */
  regraft_dist = (unsigned int *) calloc(total_edge_count, sizeof(unsigned int));
  if (!regraft_dist)
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for regraft distances\n");
    return PLL_FAILURE;
  }

  for (i = 0; i < redge_count; ++i)
    regraft_dist[i] = params->radius_min;

  regraft_edges = 0;
  j = 0;
  while ((r_edge = regraft_nodes[j]) != NULL)
  {
    /* do not re-insert back into the pruning branch */
    if (r_edge == orig_prune_edge || r_edge == orig_prune_edge->back ||
        !pllmod_treeinfo_check_constraint(treeinfo, p_edge, r_edge))
    {
      ++j;
      continue;
    }

    regraft_edges++;

    /* regraft p_edge on r_edge*/
    pllmod_treeinfo_get_branch_length_all(treeinfo, r_edge, regraft_length);

    /* distance to the current regraft edge */
    r_dist = regraft_dist[j];

    /* regraft into the candidate branch */
    retval = algo_utree_regraft(treeinfo, params, p_edge, r_edge);
    assert(retval == PLL_SUCCESS);

    /* place root at the pruning branch and invalidate CLV at the new root */
    pllmod_treeinfo_set_root(treeinfo, p_edge);
    pllmod_treeinfo_invalidate_clv(treeinfo, p_edge);

    /* save branch lengths */
    pllmod_treeinfo_get_branch_length_all(treeinfo, p_edge, b1);
    pllmod_treeinfo_get_branch_length_all(treeinfo, p_edge->next, b2);
    pllmod_treeinfo_get_branch_length_all(treeinfo, p_edge->next->next, b3);

    /* make sure branches are within limits */
    algo_unode_fix_length(treeinfo, p_edge->next, params->bl_min, params->bl_max);
    algo_unode_fix_length(treeinfo, p_edge->next->next, params->bl_min, params->bl_max);

    /* invalidate p-matrices */
    pllmod_treeinfo_invalidate_pmatrix(treeinfo, p_edge->next);
    pllmod_treeinfo_invalidate_pmatrix(treeinfo, p_edge->next->next);

    /* recompute p-matrices for branches adjacent to regrafting point */
    algo_update_pmatrix(treeinfo, p_edge->next);
    algo_update_pmatrix(treeinfo, p_edge->next->next);

    /* re-compute invalid CLVs, and get tree logLH */
    loglh = pllmod_treeinfo_compute_loglh_flex(treeinfo, 1, 0);

    if (params->thorough)
    {
      /* optimize 3 adjacent branches and get tree logLH */
      loglh = algo_optimize_bl_triplet(p_edge,
                                       treeinfo,
                                       params,
                                       1.0);

      if (!loglh)
      {
        free(regraft_nodes);
        free(regraft_dist);

        return PLL_FAILURE;
      }
    }

    if (loglh > entry->lh)
    {
//      double lh2 = pllmod_treeinfo_compute_loglh(treeinfo, 0);
//      printf("loglh part / full: %f / %f\n", loglh, lh2);

      entry->lh = loglh;
      entry->r_node = r_edge;
      pllmod_treeinfo_get_branch_length_all(treeinfo, p_edge, entry->b1);
      pllmod_treeinfo_get_branch_length_all(treeinfo, p_edge->next, entry->b2);
      pllmod_treeinfo_get_branch_length_all(treeinfo, p_edge->next->next, entry->b3);
    }

    // restore original branch lengths
    pllmod_treeinfo_set_branch_length_all(treeinfo, p_edge, b1);
    pllmod_treeinfo_set_branch_length_all(treeinfo, p_edge->next, b2);
    pllmod_treeinfo_set_branch_length_all(treeinfo, p_edge->next->next, b3);

    pllmod_treeinfo_invalidate_pmatrix(treeinfo, p_edge);
    pllmod_treeinfo_invalidate_pmatrix(treeinfo, p_edge->next);
    pllmod_treeinfo_invalidate_pmatrix(treeinfo, p_edge->next->next);

    /* rollback the REGRAFT */
    pll_unode_t * pruned_tree = pllmod_utree_prune(p_edge);
    pllmod_treeinfo_set_branch_length_all(treeinfo, pruned_tree, regraft_length);
    pllmod_treeinfo_invalidate_pmatrix(treeinfo, pruned_tree);

    /* recompute p-matrix for the pendant branch of the pruned subtree */
    algo_update_pmatrix(treeinfo, p_edge);

    /* recompute p-matrix for the "old" regraft branch */
    algo_update_pmatrix(treeinfo, pruned_tree);

    descent = r_dist < params->radius_max;
    if (cutoff_info && loglh < cutoff_info->lh_start)
    {
      cutoff_info->lh_dec_count++;
      cutoff_info->lh_dec_sum += cutoff_info->lh_start - loglh;
      descent = descent &&
                (cutoff_info->lh_start - loglh) < cutoff_info->lh_cutoff;
    }

    if (r_edge->next && descent)
    {
      regraft_nodes[redge_count] = r_edge->next->back;
      regraft_nodes[redge_count+1] = r_edge->next->next->back;
      regraft_dist[redge_count] = regraft_dist[redge_count+1] = r_dist+1;
      redge_count += 2;
    }

    ++j;
    assert(j < total_edge_count);
  }

  /* done with regrafting; restore old root */
  pllmod_treeinfo_set_root(treeinfo, orig_prune_edge);

  /* re-insert into the original pruning branch */
  retval = pllmod_utree_regraft(p_edge, orig_prune_edge);
  assert(retval == PLL_SUCCESS || (pll_errno & PLLMOD_TREE_ERROR_SPR_MASK));

  /* restore original branch length */
  pllmod_treeinfo_set_branch_length_all(treeinfo, p_edge, z1);
  pllmod_treeinfo_set_branch_length_all(treeinfo, p_edge->next, z2);
  pllmod_treeinfo_set_branch_length_all(treeinfo, p_edge->next->next, z3);

  /* invalidate p-matrices */
  pllmod_treeinfo_invalidate_pmatrix(treeinfo, p_edge);
  pllmod_treeinfo_invalidate_pmatrix(treeinfo, p_edge->next);
  pllmod_treeinfo_invalidate_pmatrix(treeinfo, p_edge->next->next);

  free(regraft_nodes);
  free(regraft_dist);

  return PLL_SUCCESS;
}

static double reinsert_nodes(pllmod_treeinfo_t * treeinfo, pll_unode_t ** nodes,
                             int node_count,
                             pllmod_rollback_list_t * rollback_list,
                             pllmod_bestnode_list_t * best_node_list,
                             cutoff_info_t * cutoff_info,
                             const pllmod_search_params_t * params)
{
  int i;

  double loglh   = pllmod_treeinfo_compute_loglh(treeinfo, 0);
  double best_lh = loglh;

  node_entry_t spr_entry;

  spr_entry.b1 = params->brlen_buf[7];
  spr_entry.b2 = params->brlen_buf[8];
  spr_entry.b3 = params->brlen_buf[9];

  pll_tree_rollback_t * rollback = rollback_list->list + rollback_list->current;

  for (i = 0; i < node_count; ++i)
  {
    pll_unode_t * p_edge = nodes[i];

    assert(!pllmod_utree_is_tip(p_edge));

    /* if remaining pruned tree would only contain 2 taxa, skip this node */
    if (pllmod_utree_is_tip(p_edge->next->back) &&
        pllmod_utree_is_tip(p_edge->next->next->back))
      continue;

    spr_entry.p_node = p_edge;

    if (cutoff_info)
      cutoff_info->lh_start = best_lh;

    int retval = best_reinsert_edge(treeinfo, &spr_entry, cutoff_info, params);
    if (!retval)
    {
      /* return and spread error */
      return 0;
    }

    pll_unode_t * best_r_edge = spr_entry.r_node;

    /* original placement is the best for the current node -> move on to the next one */
    if (!best_r_edge || best_r_edge == p_edge || best_r_edge == p_edge->back ||
        best_r_edge->back == p_edge)
    {
      continue;
    }

    /* LH improved -> re-apply the SPR move */
    if (spr_entry.lh - best_lh > 1e-6)
    {
      /* re-apply best SPR move for the node */
      DBG("SPR: %u -> (%u %u)\n", p_edge->clv_index,
          best_r_edge->clv_index, best_r_edge->back->clv_index);

      pll_unode_t * orig_prune_edge = p_edge->next->back;
      int retval = algo_utree_spr(treeinfo, params, p_edge, best_r_edge, rollback);
      assert(retval == PLL_SUCCESS);
      if (!retval)
        return PLL_FAILURE;

      algo_unode_fix_length(treeinfo, orig_prune_edge, params->bl_min, params->bl_max);

      /* increment rollback slot counter to save SPR history */
      rollback = algo_rollback_list_next(rollback_list);

      if (params->thorough)
      {
        /* restore optimized branch length */
        pllmod_treeinfo_set_branch_length_all(treeinfo, p_edge, spr_entry.b1);
        pllmod_treeinfo_set_branch_length_all(treeinfo, p_edge->next, spr_entry.b2);
        pllmod_treeinfo_set_branch_length_all(treeinfo, p_edge->next->next, spr_entry.b3);
      }
      else
      {
        /* make sure branches are within limits */
        algo_unode_fix_length(treeinfo, p_edge->next, params->bl_min, params->bl_max);
        algo_unode_fix_length(treeinfo, p_edge->next->next, params->bl_min, params->bl_max);
      }

      assert(spr_entry.lh > best_lh);

      best_lh = spr_entry.lh;

      DBG("New best: %f\n", best_lh);
    }
    else if (best_r_edge)
    {
      /* LH didn't improve but could be still high enough to be in top-20 */
      spr_entry.rollback_num = algo_rollback_list_abspos(rollback_list);

      DBG("SAVE bestnode[%d]: %lf %lf %lf\n",
          spr_entry.rollback_num, spr_entry.b1[0], spr_entry.b2[0], spr_entry.b3[0]);

      algo_bestnode_list_save(best_node_list, &spr_entry);
      loglh = spr_entry.lh;
    }

    #ifdef _ULTRACHECK
    double tmp_logh = pllmod_treeinfo_compute_loglh(treeinfo, 0);
    if(fabs(tmp_logh - best_lh) > 10e-6)
    {
      printf("%f %f\n", tmp_logh, best_lh);
      assert(0);
    }
    #endif

    DBG("LogLikelihood after SPRs for node %d (idx %d, clv %d): %f, best LH: %f\n",
         i, p_edge->node_index, p_edge->clv_index, loglh, best_lh);
  }

  return loglh;
}

PLL_EXPORT double pllmod_algo_spr_round(pllmod_treeinfo_t * treeinfo,
                                        unsigned int radius_min,
                                        unsigned int radius_max,
                                        unsigned int ntopol_keep,
                                        pll_bool_t thorough,
                                        int brlen_opt_method,
                                        double bl_min,
                                        double bl_max,
                                        int smoothings,
                                        double epsilon,
                                        cutoff_info_t * cutoff_info,
                                        double subtree_cutoff)
{
  unsigned int i;
  double loglh, best_lh;
  pllmod_search_params_t params;
  int retval;
  int brlen_unlinked;

  unsigned int allnodes_count;
  pll_unode_t ** allnodes = NULL;

  size_t rollback_slots;
  size_t toplist_slots;
  unsigned int brlen_set_count;
  pllmod_rollback_list_t * rollback_list = NULL;
  pllmod_bestnode_list_t * bestnode_list = NULL;
  pll_tree_rollback_t * rollback;
  size_t rollback_counter;
  pll_tree_rollback_t * rollback2 = NULL;
  int toplist_index;

  node_entry_t * spr_entry;
  pll_unode_t * p_edge, * r_edge;

  pllmod_treeinfo_topology_t * best_topol = NULL;
#ifndef  PLLMOD_SEARCH_GREEDY_BLO
  pllmod_treeinfo_topology_t * tmp_topol = NULL;
#endif

  double static_brlen_buf[BRLEN_BUF_COUNT];

  /* process search params */
  params.thorough = thorough;
  params.ntopol_keep = ntopol_keep;
  params.radius_min = radius_min;
  params.radius_max = radius_max;
  params.bl_min = bl_min;
  params.bl_max = bl_max;
  params.smoothings = smoothings;
  params.brlen_opt_method = brlen_opt_method;

  brlen_unlinked = (treeinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED) ? 1 : 0;

  /* reset error */
  pll_errno = 0;

  /* allocate brlen buffers */
  for (i = 0; i < BRLEN_BUF_COUNT; ++i)
  {
    params.brlen_buf[i] = brlen_unlinked ?
            (double *) calloc(treeinfo->init_partition_count, sizeof(double)) :
            &static_brlen_buf[i];
  }

  /* allocate rollback_info slots */
  rollback_slots = params.ntopol_keep;
  rollback_list = algo_rollback_list_create(rollback_slots);
  if (!rollback_list)
  {
    /* return and spread error */
    goto error_exit;
  }

  /* allocate best node slots */
  toplist_slots = params.thorough ? params.ntopol_keep : params.ntopol_keep * 3;
  brlen_set_count = brlen_unlinked ? treeinfo->init_partition_count : 1;
  bestnode_list = algo_bestnode_list_create(toplist_slots, brlen_set_count);
  if (!bestnode_list)
  {
    /* return and spread error */
    goto error_exit;
  }

  if (cutoff_info)
  {
    cutoff_info->lh_dec_count = 0;
    cutoff_info->lh_dec_sum = 0.;
  }

  loglh   = pllmod_treeinfo_compute_loglh(treeinfo, 0);
  best_lh = loglh;

  /* query all nodes */
  allnodes_count = (treeinfo->tip_count - 2) * 3;
  allnodes = (pll_unode_t **) calloc (allnodes_count, sizeof(pll_unode_t *));
  if (!allnodes)
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory nodes list\n");
    goto error_exit;
  }

  unsigned int node_count = algo_query_allnodes(treeinfo->root, allnodes);
  assert(node_count == allnodes_count);

  loglh = reinsert_nodes(treeinfo,
                         allnodes,
                         allnodes_count,
                         rollback_list,
                         bestnode_list,
                         cutoff_info,
                         &params);
  if (!loglh)
  {
    /* return and spread error */
    goto error_exit;
  }

  /* in FAST mode, we re-insert a subset of best-scoring subtrees with BLO
   * (i.e., in SLOW mode) */
  if (!params.thorough && bestnode_list->current > 0)
  {
    params.thorough = PLL_TRUE;
    for (i = 0; bestnode_list->list[i].p_node != NULL; i++)
    {
      allnodes[i] = bestnode_list->list[i].p_node;
      bestnode_list->list[i].p_node = NULL;
    }

    DBG("\nThorough re-insertion of %u best-scoring nodes...\n", i);

    loglh = reinsert_nodes(treeinfo,
                           allnodes,
                           i,
                           rollback_list,
                           bestnode_list,
                           cutoff_info,
                           &params);
    if (!loglh)
    {
      /* return and spread error */
      goto error_exit;
    }
  }

  free(allnodes);
  allnodes = NULL;

  best_lh = algo_optimize_bl_all(treeinfo,
                                 &params,
                                 epsilon,
                                 0.25);
  DBG("Best tree LH after BLO: %f\n", best_lh);

  if (!best_lh)
  {
    /* return and spread error */
    goto error_exit;
  }

  best_topol = pllmod_treeinfo_get_topology(treeinfo, NULL);
  if (!best_topol)
    goto error_exit;

  /* Restore best topologies and re-evaluate them after full BLO.
  NOTE: some SPRs were applied (if they improved LH) and others weren't.
  Therefore in order to restore the original topology, we need to either rollback
  an SPR (if it was already applied), or re-do it again (if it wasn't applied).
  We perform it by simultaneously iterating over the history of applied SPRs
  (rollback_list) and over the list of not-applied SPRs which resulted
  in topologies with the highest LH (bestnode_list).
  */
  rollback_counter = 0;
  toplist_index = -1;
  rollback2 = (pll_tree_rollback_t *) calloc(1, sizeof(pll_tree_rollback_t));
  if (!rollback2)
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for additional rollback list\n");
    goto error_exit;
  }
  int undo_SPR = 0;

#ifdef DEBUG
  algo_bestnode_list_print(bestnode_list);
#endif

  while (rollback_counter < rollback_list->size)
  {
    const int rollback_num = algo_rollback_list_abspos(rollback_list);
    toplist_index = algo_bestnode_list_next_index(bestnode_list,
                                                  rollback_num,
                                                  toplist_index);

    if (toplist_index == -1)
    {
      /* no more topologies for this rollback, so we go one slot back */
      rollback = algo_rollback_list_prev(rollback_list);

      if (!rollback || !rollback->SPR.prune_edge)
      {
        DBG("  Rollback slot %d is empty, exiting the loop...\n",
            rollback_list->current);
        break;
      }

      DBG("  Rollback BL: %.12lf %.12lf %.12lf %.12lf\n\n", rollback->SPR.prune_bl,
          rollback->SPR.prune_left_bl, rollback->SPR.prune_right_bl, rollback->SPR.regraft_bl);

      DBG("  Undoing SPR %lu (slot %d)... ", rollback_counter,
          rollback_list->current);

      retval = pllmod_tree_rollback(rollback);
      assert(retval == PLL_SUCCESS);

      rollback_counter++;

      undo_SPR = 0;
    }
    else
    {
      if ((unsigned int) toplist_index > params.ntopol_keep)
        continue;

      spr_entry = &bestnode_list->list[toplist_index];
      p_edge = spr_entry->p_node;
      r_edge = spr_entry->r_node;

      if (!p_edge)
      {
        DBG("    SPR slot %d is empty, exiting the loop...\n", toplist_index);
        break;
      }
      else
      {
        DBG("    Evaluating topology %d (idx %d, clv %d -> idx %d, clv %d), old LH: %f... ",
            toplist_index,
            p_edge->node_index,
            p_edge->clv_index,
            r_edge->node_index,
            r_edge->clv_index,
            spr_entry->lh);
      }

      if (!pllmod_treeinfo_check_constraint(treeinfo, p_edge, r_edge))
      {
        DBG("Topological constraint check failed, skip the topology.\n");
        continue;
      }

      /* re-apply best SPR move for the node */
      retval = pllmod_utree_spr(p_edge, r_edge, rollback2);
      assert(retval == PLL_SUCCESS);

#ifndef  PLLMOD_SEARCH_GREEDY_BLO
      /* save topology with original branch length before BLO */
      tmp_topol = pllmod_treeinfo_get_topology(treeinfo, tmp_topol);
#endif

      /* make sure that original prune branch length does not exceed maximum */
      algo_unode_fix_length(treeinfo, rollback2->SPR.regraft_edge, params.bl_min, params.bl_max);

      /* restore optimized branch lengths */
      if (thorough)
      {
        pllmod_treeinfo_set_branch_length_all(treeinfo, p_edge, spr_entry->b1);
        pllmod_treeinfo_set_branch_length_all(treeinfo, p_edge->next, spr_entry->b2);
        pllmod_treeinfo_set_branch_length_all(treeinfo, p_edge->next->next, spr_entry->b3);

        DBG("RESTORE bl triplet[%d/%d]: %lf %lf %lf\n", toplist_index,
            bestnode_list->brlen_set_count, spr_entry->b1[0], spr_entry->b2[0], spr_entry->b3[0]);
      }
      else
      {
        algo_unode_fix_length(treeinfo, p_edge, params.bl_min, params.bl_max);
        algo_unode_fix_length(treeinfo, p_edge->next, params.bl_min, params.bl_max);
        algo_unode_fix_length(treeinfo, p_edge->next->next, params.bl_min, params.bl_max);
      }

      undo_SPR = 1;
    }

    /* now optimize all the branches */
    double loglh;
    loglh = algo_optimize_bl_all(treeinfo,
                                 &params,
                                 epsilon,
                                 0.25);

    if (!loglh)
    {
      /* return and spread error */
      goto error_exit;
    }

    DBG("  new LH after BLO: %f\n", loglh);
    assert(loglh > -INFINITY);

    if (loglh - best_lh > 0.01)
    {
      DBG("Best tree LH: %f\n", loglh);

      best_topol = pllmod_treeinfo_get_topology(treeinfo, best_topol);
      if (!best_topol)
        goto error_exit;

      best_lh = loglh;
    }

    if (undo_SPR)
    {
#ifndef  PLLMOD_SEARCH_GREEDY_BLO
      /* restore original brlens */
      retval = pllmod_treeinfo_set_topology(treeinfo, tmp_topol);
      if (!retval)
        goto error_exit;
#endif

      /* rollback the SPR */
      retval = pllmod_tree_rollback(rollback2);
      assert(retval == PLL_SUCCESS);
    }
  }

  free(rollback2);

  algo_bestnode_list_destroy(bestnode_list);
  algo_rollback_list_destroy(rollback_list);

  if (brlen_unlinked)
  {
    for (i = 0; i < BRLEN_BUF_COUNT; ++i)
      free(params.brlen_buf[i]);
  }

  /* update LH cutoff */
  if (cutoff_info)
  {
    cutoff_info->lh_cutoff =
        subtree_cutoff * (cutoff_info->lh_dec_sum / cutoff_info->lh_dec_count);
  }

  if (best_topol)
  {
    retval = pllmod_treeinfo_set_topology(treeinfo, best_topol);
    pllmod_treeinfo_destroy_topology(best_topol);
    if (!retval)
      goto error_exit;
  }

#ifndef  PLLMOD_SEARCH_GREEDY_BLO
  if (tmp_topol)
    pllmod_treeinfo_destroy_topology(tmp_topol);
#endif

  /* update partials and CLVs */
  loglh = pllmod_treeinfo_compute_loglh(treeinfo, 0);
  if (fabs(loglh - best_lh) > 1e-6)
  {
    printf("LH mismatch: %.12f  != %.12f\n", best_lh, loglh);
    assert(fabs(loglh - best_lh) < 1e-6);
  }

  return loglh;

error_exit:
  /* cleanup */
  if (allnodes)
    free(allnodes);
  if (rollback2)
    free(rollback2);
  algo_bestnode_list_destroy(bestnode_list);
  algo_rollback_list_destroy(rollback_list);

  /* make sure libpll error code is set and exit */
  assert(pll_errno);
  return 0;
}
