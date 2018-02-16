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

typedef struct spr_params
{
  int thorough;
  int radius_min;
  int radius_max;
  int ntopol_keep;
  double bl_min;
  double bl_max;
  int smoothings;
  int brlen_opt_method;
} pllmod_search_params_t;

typedef struct rollback_list
{
  pll_tree_rollback_t * list;
  int current;
  size_t round;
  size_t size;
} pllmod_rollback_list_t;

typedef struct node_entry {
  pll_unode_t * p_node;
  pll_unode_t * r_node;
  double b1, b2, b3;
  double lh;
  unsigned int rollback_num;
} node_entry_t;

typedef struct bestnode_list
{
  node_entry_t * list;
  int current;
  size_t size;
} pllmod_bestnode_list_t;

static void algo_query_allnodes_recursive(pll_unode_t * node,
                                          pll_unode_t ** buffer,
                                          int * index)
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

static int algo_query_allnodes(pll_unode_t * root, pll_unode_t ** buffer)
{
  assert(root && buffer);

  int index = 0;

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
  else if (rollback_list->round > 0)
  {
    rollback_list->round--;
    rollback_list->current = rollback_list->size - 1;
  }
  else
    return NULL;

  return rollback_list->list + rollback_list->current;
}

static pll_tree_rollback_t * algo_rollback_list_next(
                               pllmod_rollback_list_t * rollback_list)
{
  if (rollback_list->current < ((int) rollback_list->size - 1))
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

static pllmod_bestnode_list_t * algo_bestnode_list_create(size_t slots)
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
  }
  return bestnode_list;
}

static void algo_bestnode_list_destroy(pllmod_bestnode_list_t * bestnode_list)
{
 if (bestnode_list)
 {
   if (bestnode_list->list)
     free(bestnode_list->list);
   free(bestnode_list);
 }
}

static void algo_bestnode_list_copy_entry(node_entry_t * dst,
                                          const node_entry_t * src)
{
  memcpy(dst, src, sizeof(node_entry_t));
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
    algo_bestnode_list_copy_entry(&list[j], &list[j-1]);
    --j;
  }

  assert(idx >= 0 && idx < list_size);

  algo_bestnode_list_copy_entry(&list[idx], entry);
}

static int algo_bestnode_list_next_index(pllmod_bestnode_list_t * best_node_list,
                                         unsigned int rollback_num,
                                         int curr_index)
{
  assert(curr_index >= -1);

  node_entry_t * list = best_node_list->list;
  const int list_size = best_node_list->size;

  do
  {
    curr_index++;
    if (curr_index >= (int) list_size || !list[curr_index].p_node)
      return -1;
  }
  while (list[curr_index].rollback_num != rollback_num);

  return curr_index;
}

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

static void algo_unode_fix_length(pll_unode_t * node, double bl_min, double bl_max)
{
  if (node->length < bl_min)
    pllmod_utree_set_length(node, bl_min);
  else if (node->length > bl_max)
    pllmod_utree_set_length(node, bl_max);
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

      // TODO: extend for unlinked per-partition branch lengths
      double p_brlen = edge->length;
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
  double regraft_length;
  int r_dist;
  double z1, z2, z3;
  unsigned int redge_count = 0;
  unsigned int ncount;
  int retval;
  int * regraft_dist;
  int descent;
  double b1, b2, b3;
  double loglh;

  pll_unode_t * p_edge = entry->p_node;
  const size_t total_edge_count = 2 * treeinfo->tip_count - 3;

  entry->r_node = NULL;
  entry->lh = PLLMOD_OPT_LNL_UNLIKELY;

  /* save original branch lengths at the pruning point */
  z1 = p_edge->length;
  z2 = p_edge->next->length;
  z3 = p_edge->next->next->length;

  pllmod_treeinfo_set_root(treeinfo, p_edge);

  /* recompute all CLVs and p-matrices before pruning */
  pllmod_treeinfo_compute_loglh(treeinfo, 0);

  orig_prune_edge = pllmod_utree_prune(p_edge);
  if (!orig_prune_edge)
  {
    /* check that errno was set correctly */
    assert(pll_errno & PLLMOD_TREE_ERROR_SPR_MASK);
    return PLL_FAILURE;
  }

  algo_unode_fix_length(orig_prune_edge, params->bl_min, params->bl_max);

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
  regraft_dist = (int *) calloc(total_edge_count, sizeof(int));
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
    if (r_edge == orig_prune_edge || r_edge == orig_prune_edge->back)
      continue;

    regraft_edges++;

    /* regraft p_edge on r_edge*/
    regraft_length = r_edge->length;

    /* distance to the current regraft edge */
    r_dist = regraft_dist[j];

    /* regraft into the candidate branch */
    retval = pllmod_utree_regraft(p_edge, r_edge);
    assert(retval == PLL_SUCCESS);

    /* invalidate p-matrices */
    pllmod_treeinfo_invalidate_pmatrix(treeinfo, p_edge->next);
    pllmod_treeinfo_invalidate_pmatrix(treeinfo, p_edge->next->next);

    /* place root at the pruning branch and invalidate CLV at the new root */
    pllmod_treeinfo_set_root(treeinfo, p_edge);
    pllmod_treeinfo_invalidate_clv(treeinfo, p_edge);

    /* save branch lengths */
    b1 = p_edge->length;
    b2 = p_edge->next->length;
    b3 = p_edge->next->next->length;

    /* make sure branches are within limits */
    algo_unode_fix_length(p_edge->next, params->bl_min, params->bl_max);
    algo_unode_fix_length(p_edge->next->next, params->bl_min, params->bl_max);

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
      entry->lh = loglh;
      entry->r_node = r_edge;
      entry->b1 = p_edge->length;
      entry->b2 = p_edge->next->length;
      entry->b3 = p_edge->next->next->length;
    }

    // restore original branch lengths
    pllmod_utree_set_length(p_edge, b1);
    pllmod_utree_set_length(p_edge->next, b2);
    pllmod_utree_set_length(p_edge->next->next, b3);

    pllmod_treeinfo_invalidate_pmatrix(treeinfo, p_edge);
    pllmod_treeinfo_invalidate_pmatrix(treeinfo, p_edge->next);
    pllmod_treeinfo_invalidate_pmatrix(treeinfo, p_edge->next->next);

    /* rollback the REGRAFT */
    pll_unode_t * pruned_tree = pllmod_utree_prune(p_edge);
    pllmod_utree_set_length(pruned_tree, regraft_length);
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
  pllmod_utree_set_length(p_edge, z1);
  pllmod_utree_set_length(p_edge->next, z2);
  pllmod_utree_set_length(p_edge->next->next, z3);

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
      pll_unode_t * orig_prune_edge = p_edge->next->back;
      int retval = pllmod_utree_spr(p_edge, best_r_edge, rollback);
      assert(retval == PLL_SUCCESS);

      algo_unode_fix_length(orig_prune_edge, params->bl_min, params->bl_max);

      /* increment rollback slot counter to save SPR history */
      rollback = algo_rollback_list_next(rollback_list);

      if (params->thorough)
      {
        /* restore optimized branch length */
        pllmod_utree_set_length(p_edge, spr_entry.b1);
        pllmod_utree_set_length(p_edge->next, spr_entry.b2);
        pllmod_utree_set_length(p_edge->next->next, spr_entry.b3);
      }
      else
      {
        /* make sure branches are within limits */
        algo_unode_fix_length(p_edge->next, params->bl_min, params->bl_max);
        algo_unode_fix_length(p_edge->next->next, params->bl_min, params->bl_max);
      }

      assert(spr_entry.lh > best_lh);

      best_lh = spr_entry.lh;

      DBG("New best: %f\n", best_lh);
    }
    else if (best_r_edge)
    {
      /* LH didn't improve but could be still high enough to be in top-20 */
      spr_entry.rollback_num = algo_rollback_list_abspos(rollback_list);
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
                                        int radius_min,
                                        int radius_max,
                                        int ntopol_keep,
                                        int thorough,
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

  int allnodes_count;
  pll_unode_t ** allnodes = NULL;

  size_t rollback_slots;
  size_t toplist_slots;
  pllmod_rollback_list_t * rollback_list = NULL;
  pllmod_bestnode_list_t * bestnode_list = NULL;
  pll_tree_rollback_t * rollback;
  size_t rollback_counter;
  pll_tree_rollback_t * rollback2 = NULL;
  int toplist_index;

  node_entry_t * spr_entry;
  pll_unode_t * p_edge, * r_edge;

  pll_unode_t * best_tree;
#ifndef  PLLMOD_SEARCH_GREEDY_BLO
  pll_unode_t * tmp_tree = NULL;
#endif

  /* process search params */
  params.thorough = thorough;
  params.ntopol_keep = ntopol_keep;
  params.radius_min = radius_min;
  params.radius_max = radius_max;
  params.bl_min = bl_min;
  params.bl_max = bl_max;
  params.smoothings = smoothings;
  params.brlen_opt_method = brlen_opt_method;

  /* reset error */
  pll_errno = 0;

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
  bestnode_list = algo_bestnode_list_create(toplist_slots);
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
  allnodes = (pll_unode_t **) calloc ((size_t) allnodes_count,
                                      sizeof(pll_unode_t *));
  if (!allnodes)
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory nodes list\n");
    goto error_exit;
  }

  assert(algo_query_allnodes(treeinfo->root, allnodes) == allnodes_count);

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
    params.thorough = 1;
    for (i = 0; bestnode_list->list[i].p_node != NULL; i++)
    {
      allnodes[i] = bestnode_list->list[i].p_node;
      bestnode_list->list[i].p_node = NULL;
    }

    DBG("\nThorough re-insertion of %d best-scoring nodes...\n", i);

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

  best_tree = pll_utree_graph_clone(treeinfo->root);

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
      if (toplist_index > params.ntopol_keep)
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

      /* re-apply best SPR move for the node */
      retval = pllmod_utree_spr(p_edge, r_edge, rollback2);
      assert(retval == PLL_SUCCESS);

#ifndef  PLLMOD_SEARCH_GREEDY_BLO
      /* clone the tree before BLO to save original branch length */
      // TODO: make it more efficient
      tmp_tree = treeinfo->root;
      treeinfo->root = pll_utree_graph_clone(treeinfo->root);
#endif

      /* make sure that original prune branch length does not exceed maximum */
      algo_unode_fix_length(rollback2->SPR.regraft_edge, params.bl_min, params.bl_max);

      /* restore optimized branch lengths */
      pllmod_utree_set_length(p_edge, spr_entry->b1);
      pllmod_utree_set_length(p_edge->next, spr_entry->b2);
      pllmod_utree_set_length(p_edge->next->next, spr_entry->b3);

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

      if (best_tree)
        pll_utree_graph_destroy(best_tree, NULL);

      best_tree = pll_utree_graph_clone(treeinfo->root);
      best_lh = loglh;
    }

    if (undo_SPR)
    {
#ifndef  PLLMOD_SEARCH_GREEDY_BLO
      /* restore original brlens */
      pll_utree_graph_destroy(treeinfo->root, NULL);
      treeinfo->root = tmp_tree;
#endif

      /* rollback the SPR */
      retval = pllmod_tree_rollback(rollback2);
      assert(retval == PLL_SUCCESS);
    }
  }

  if (best_tree)
  {
    pll_utree_graph_destroy(treeinfo->root, NULL);
    treeinfo->root = best_tree;
  }

  free(rollback2);

  algo_bestnode_list_destroy(bestnode_list);
  algo_rollback_list_destroy(rollback_list);

  /* update LH cutoff */
  if (cutoff_info)
  {
    cutoff_info->lh_cutoff =
        subtree_cutoff * (cutoff_info->lh_dec_sum / cutoff_info->lh_dec_count);
  }

  /* update partials and CLVs */
  loglh = pllmod_treeinfo_compute_loglh(treeinfo, 0);
  assert(fabs(loglh - best_lh) < 1e-6);

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
