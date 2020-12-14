/*
 Copyright (C) 2016 Diego Darriba

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

 Contact: Diego Darriba <Diego.Darriba@h-its.org>,
 Exelixis Lab, Heidelberg Instutute for Theoretical Studies
 Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
 */

 /**
  * @file pll_tree.c
  *
  * @brief Operations on tree structures
  *
  * @author Diego Darriba
  */

#include "pll_tree.h"
#include "tree_hashtable.h"
#include "../pllmod_common.h"

#define UNIMPLEMENTED 0

static int rtree_rollback_tbr(pll_tree_rollback_t * rollback_info);
static int rtree_rollback_spr(pll_tree_rollback_t * rollback_info);
static int rtree_rollback_nni(pll_tree_rollback_t * rollback_info);
static int utree_rollback_tbr(pll_tree_rollback_t * rollback_info);
static int utree_rollback_spr(pll_tree_rollback_t * rollback_info);
static int utree_rollback_nni(pll_tree_rollback_t * rollback_info);
static int utree_find_node_in_subtree(pll_unode_t * root, pll_unode_t * node);
static int cb_update_matrices_partials(pll_unode_t * node, void *data);
static void shuffle_tree_nodes(const pll_utree_t * tree, unsigned int seed);
static void split_multi_node(pll_utree_t * tree, pll_unode_t * first,
                             pll_unode_t * last, unsigned int degree);
static char * default_support_fmt(double support);

struct cb_params
{
  const unsigned int * params_indices;
  pll_partition_t * partition;
  int update_pmatrices;
  int update_partials;
};

/******************************************************************************/
/* Topological rearrangements */

/**
 * Performs one TBR move by applying a bisection and a reconnection.
 * The CLV, scaler and pmatrix indices are updated.
 *
 * @param[in] b_edge bisection point
 * @param[in] r_edge reconnection point
 * @param[out] rollback_info Rollback information for undoing this move.
 *                           If it is NULL, rollback information is ignored.
 *
 * @return PLL_SUCCESS if the move was applied correctly,
 *         PLL_FAILURE otherwise (check pll_errmsg for details)
 */
PLL_EXPORT int pllmod_utree_tbr(pll_unode_t * b_edge,
                                pll_tree_edge_t * r_edge,
                                pll_tree_rollback_t * rollback_info)
{
  pll_unode_t *parent, *child;

  /* validate if the move can be applied */

  /* 1. bisection point must not be a leaf branch */
  if (!(b_edge->next && b_edge->back->next))
  {
    pllmod_set_error(PLLMOD_TREE_ERROR_TBR_LEAF_BISECTION,
                     "attempting to bisect at a leaf node");
    return PLL_FAILURE;
  }

  /* 2. reconnection edges are different from bisection point */
  if (b_edge == r_edge->edge.utree.parent ||
      b_edge == r_edge->edge.utree.parent->back ||
      b_edge == r_edge->edge.utree.child ||
      b_edge == r_edge->edge.utree.child->back ||
      b_edge->back == r_edge->edge.utree.parent ||
      b_edge->back == r_edge->edge.utree.parent->back ||
      b_edge->back == r_edge->edge.utree.child ||
      b_edge->back == r_edge->edge.utree.child->back)
  {
    pllmod_set_error(PLLMOD_TREE_ERROR_TBR_OVERLAPPED_NODES,
                     "TBR nodes are overlapped");
    return PLL_FAILURE;
  }

  /* 3. reconnection edges must belong to different subtrees rooted at b_edge
   *    and b_edge->back
   */
  if (!(utree_find_node_in_subtree(b_edge, r_edge->edge.utree.parent) &&
        utree_find_node_in_subtree(b_edge->back, r_edge->edge.utree.child)) &&
      !(utree_find_node_in_subtree(b_edge->back, r_edge->edge.utree.parent) &&
        utree_find_node_in_subtree(b_edge, r_edge->edge.utree.child)))
  {
    pllmod_set_error(PLLMOD_TREE_ERROR_TBR_SAME_SUBTREE,
                     "TBR reconnection in same subtree");
    return PLL_FAILURE;
  }

  /* save rollback information */
  if (rollback_info)
  {
    rollback_info->rearrange_type     = PLLMOD_TREE_REARRANGE_TBR;
    rollback_info->rooted             = 0;
    rollback_info->TBR.bisect_edge    = (void *) b_edge;
    rollback_info->TBR.reconn_edge.edge.utree.parent = b_edge->next->next;
    rollback_info->TBR.reconn_edge.edge.utree.child  = b_edge->back->next->next;
    rollback_info->TBR.reconn_edge.length = b_edge->length;

    rollback_info->TBR.bisect_left_bl = r_edge->edge.utree.parent->length;
    rollback_info->TBR.bisect_right_bl = r_edge->edge.utree.child->length;

    rollback_info->TBR.reconn_parent_left_bl  = b_edge->next->length;
    rollback_info->TBR.reconn_parent_right_bl = b_edge->next->next->length;
    rollback_info->TBR.reconn_child_left_bl   = b_edge->back->next->length;
    rollback_info->TBR.reconn_child_right_bl  = b_edge->back->next->next->length;
  }

  /* bisect at b_edge */
  pllmod_utree_bisect(b_edge, &parent, &child);

  /* reconnect at r_edge */
  pllmod_utree_reconnect(r_edge,
                      b_edge);

  return PLL_SUCCESS;
}



/**
 * Performs one SPR move
 * The CLV, scaler and pmatrix indices are updated.
 *
 * @param[in] p_edge Edge to be pruned
 * @param[in] r_edge Edge to be regrafted
 * @param[out] rollback_info Rollback information for undoing this move.
 *                           If it is NULL, rollback information is ignored.
 *
 * @return PLL_SUCCESS if the move was applied correctly,
 *         PLL_FAILURE otherwise (check pll_errmsg for details)
 */
PLL_EXPORT int pllmod_utree_spr(pll_unode_t * p_edge,
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

  retval = pll_utree_spr(p_edge,
                         r_edge,
                         0, 0, 0);

  return retval;
}

/**
 * Performs one NNI move
 * The CLV, scaler and pmatrix indices are updated.
 *
 * @param[in] edge NNI interchange edge
 * @param[in] type move type: PLL_NNI_LEFT, PLL_NNI_RIGHT
 * @param[out] rollback_info Rollback information for undoing this move.
 *                           If it is NULL, rollback information is ignored.
 *
 * @return PLL_SUCCESS if the move was applied correctly,
 *         PLL_FAILURE otherwise (check pll_errmsg for details)
 */
PLL_EXPORT int pllmod_utree_nni(pll_unode_t * edge,
                                int type,
                                pll_tree_rollback_t * rollback_info)
{
  /* validate preconditions */
  assert(edge && edge->back);

  if (!(type == PLL_UTREE_MOVE_NNI_LEFT || type == PLL_UTREE_MOVE_NNI_RIGHT))
  {
    /* invalid move */
    pllmod_set_error(PLLMOD_TREE_ERROR_NNI_INVALID_MOVE,
                     "Invalid NNI move type");
    return PLL_FAILURE;
  }
  if (pllmod_utree_is_tip(edge) || pllmod_utree_is_tip(edge->back))
  {
    /* invalid move */
    pllmod_set_error(PLLMOD_TREE_ERROR_INTERCHANGE_LEAF,
                     "Attempting to apply NNI on a leaf branch");
    return PLL_FAILURE;
  }

  /* save rollback information */
  if (rollback_info)
  {
    rollback_info->rearrange_type     = PLLMOD_TREE_REARRANGE_NNI;
    rollback_info->rooted             = 0;
    rollback_info->NNI.edge           = (void *) edge;
    rollback_info->NNI.type           = type;
    rollback_info->NNI.left_left_bl   = edge->next->length;
    rollback_info->NNI.left_right_bl  = edge->next->next->length;
    rollback_info->NNI.right_left_bl  = edge->back->next->length;
    rollback_info->NNI.right_right_bl = edge->back->next->next->length;
    rollback_info->NNI.edge_bl        = edge->length;
  }

  if (!pll_utree_nni(edge, type, 0))
    return PLL_FAILURE;

  return PLL_SUCCESS;
}

/**
 * Rollback the previous move
 * @param  rollback_info the rollback info returned by the previous move
 * @return PLL_SUCCESS if the rollback move was applied correctly,
 *         PLL_FAILURE otherwise (check pll_errmsg for details)
 */
PLL_EXPORT int pllmod_tree_rollback(pll_tree_rollback_t * rollback_info)
{
  int retval = PLL_FAILURE;
  switch (rollback_info->rearrange_type)
  {
    case PLLMOD_TREE_REARRANGE_TBR:
      {
        if (rollback_info->rooted)
          retval = rtree_rollback_tbr (rollback_info);
        else
          retval = utree_rollback_tbr (rollback_info);
      }
      break;
    case PLLMOD_TREE_REARRANGE_SPR:
      {
        if (rollback_info->rooted)
          retval = rtree_rollback_spr (rollback_info);
        else
          retval = utree_rollback_spr (rollback_info);
      }
      break;
    case PLLMOD_TREE_REARRANGE_NNI:
      {
        if (rollback_info->rooted)
          retval = rtree_rollback_nni (rollback_info);
        else
          retval = utree_rollback_nni (rollback_info);
      }
      break;
    default:
      /* unimplemented */
      assert(0);
      break;
  }
  return retval;
}



/******************************************************************************/
/* Tree construction */

PLL_EXPORT pll_utree_t * pllmod_utree_resolve_multi(const pll_utree_t * multi_tree,
                                                    unsigned int random_seed,
                                                    int * clv_index_map)
{
  if (!multi_tree)
  {
    pllmod_set_error(PLL_ERROR_PARAM_INVALID, "Parameter multi_tree is NULL.");
    return NULL;
  }

  if (multi_tree->vroot->next &&
      multi_tree->vroot->next->next == multi_tree->vroot)
  {
    pllmod_set_error(PLLMOD_TREE_ERROR_INVALID_TREE,
                     "Unrooted tree is expected but a rooted tree was provided.");
    return NULL;
  }

  pll_utree_t * bin_tree = pll_utree_clone(multi_tree);

  unsigned int tip_count = bin_tree->tip_count;
  unsigned int multi_node_count = bin_tree->tip_count + bin_tree->inner_count;
  unsigned int bin_node_count = 2 * tip_count -2;

  // 1:1 CLV index mapping for existing nodes
  if (clv_index_map)
  {
    for (unsigned int i = 0; i < multi_node_count; ++i)
    {
      const unsigned int clv_id = bin_tree->nodes[i]->clv_index;
      clv_index_map[clv_id] = (int) clv_id;
    }
  }

  if (bin_tree->binary)
    return bin_tree;

  if (random_seed)
    shuffle_tree_nodes(bin_tree, random_seed);

  bin_tree->nodes = (pll_unode_t **)realloc(bin_tree->nodes,
                                            bin_node_count*sizeof(pll_unode_t *));

  // iterate over inner nodes, resolve multifurcations, and map new->old CLV indices
  unsigned int old_inner_count = bin_tree->inner_count;
  for (unsigned int i = tip_count; i < multi_node_count; ++i)
  {
    pll_unode_t * start = bin_tree->nodes[i];
    pll_unode_t * end = NULL;
    pll_unode_t * snode = start;
    unsigned int degree = 0;
    do
    {
      end = snode;
      snode = snode->next;
      degree++;
    }
    while (snode && snode != start);

    split_multi_node(bin_tree, start, end, degree);

    assert(bin_tree->inner_count == old_inner_count + degree-3);
    if (clv_index_map)
    {
      for (unsigned int j = old_inner_count; j < bin_tree->inner_count; ++j)
      {
        unsigned int new_clv_id = bin_tree->nodes[tip_count + j]->clv_index;
        clv_index_map[new_clv_id] = (int) start->clv_index;
      }
    }
    old_inner_count = bin_tree->inner_count;
  }

  assert(bin_tree->inner_count == bin_tree->tip_count - 2);

  bin_tree->binary = 1;

  /* re-assign node indices such that:
   * (1) all 3 pll_unode's of an inner node have consecutive indices: (x, x+1, x+2)
   * (2) for any two random multifurcation resolutions R1 and R2 holds
   *      (x, x+1, x+2) in R1 iff (x, x+1, x+2) in R2
   */
  unsigned int max_node_index = tip_count;
  for (unsigned int i = tip_count; i < bin_node_count; ++i)
  {
    pll_unode_t * node = bin_tree->nodes[i];
    node->node_index = max_node_index++;
    node->next->node_index = max_node_index++;
    node->next->next->node_index = max_node_index++;
  }
  assert(max_node_index == bin_tree->tip_count + 3*bin_tree->inner_count);

  return bin_tree;
}

static pll_unode_t * unode_prev(pll_unode_t * node)
{
  if (node->next)
  {
    pll_unode_t * prev = node;
    while (prev->next != node) prev = prev->next;
    return prev;
  }
  else
    return NULL;
}

/* This function removes a branch between lnode and lnode->back by
 * "dissolving" a roundabout ("inner node triplet") that contains lnode->back,
 * and merging its remainders into the "left" roundabout as show below:
 *
 *     *-l2-*          *-r2-*                  *-l2-*
      /      \        /      \                /      \
 * --l1  x1  l3------r1  x2  r3--   ---->  --l1  x1  r2--
 *    \      /        \      /                \      /
 *     *----*          *----*                  *-r3-*
 *
 *  where l3 = lnode, r1 = lnode->back
 */
static int remove_branch(pll_unode_t * lnode)
{
  pll_unode_t * rnode = lnode->back;
  
  /* can only remove a branch between two inner nodes */
  if (!lnode->next || !rnode->next)
    return PLL_FAILURE;

  pll_unode_t * lnode_prev = unode_prev(lnode);
  pll_unode_t * lnode_next = lnode->next;
  pll_unode_t * rnode_prev = unode_prev(rnode);
  pll_unode_t * rnode_next = rnode->next;

  /* merge remaining subnodes of left and right nodes */
  lnode_prev->next = rnode_next;
  rnode_prev->next = lnode_next;

  /* update clv_index and scaler_index in right node remainder */
  while (rnode_next != lnode_next)
  {
    rnode_next->clv_index = lnode_prev->clv_index;
    rnode_next->scaler_index = lnode_prev->scaler_index;
    rnode_next->label = lnode_prev->label;
    rnode_next = rnode_next->next;
  }

  /* destroy both subnodes adjacent to the deleted branch */
  free(lnode);
  free(rnode->label);
  free(rnode);

  return PLL_SUCCESS;
}

PLL_EXPORT int pllmod_utree_collapse_branches(pll_utree_t * tree,
                                              double min_brlen)
{
  if (!tree || !tree->vroot)
  {
    pllmod_set_error(PLL_ERROR_PARAM_INVALID, "Empty tree specified!");
    return PLL_FAILURE;
  }

  double brlen_cutoff = min_brlen + PLL_ONE_EPSILON;
  unsigned int tip_count = tree->tip_count;
  unsigned int inner_count = tree->inner_count;
  unsigned int node_count = inner_count + tip_count;
  unsigned int removed_count = 0;
  unsigned int * clv2pos_map = (unsigned int *) calloc(inner_count, sizeof(unsigned int));

  if (!clv2pos_map)
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for clv2pos map!");
    return PLL_FAILURE;
  }

  /* to avoid making assumptions about node ordering in tree->nodes,
   * we build this map indexed based on clv_index */
  for (unsigned int i = tip_count; i < node_count; ++i)
  {
    unsigned int inner_clv_idx = tree->nodes[i]->clv_index - tip_count;
    clv2pos_map[inner_clv_idx] = i;
  }

  for (unsigned int i = tip_count; i < node_count; ++i)
  {
    pll_unode_t * node  = tree->nodes[i];

    /* this node has been removed in a previous iteration -> skip */
    if (!node)
      continue;

    assert (!pllmod_utree_is_tip(node));

    pll_unode_t * start_node = NULL;
    do
    {
      pll_unode_t * anode = node->back;
      if (pllmod_utree_is_tip(anode) || node->length > brlen_cutoff)
      {
        if (!start_node)
          start_node = node;
        node = node->next;
      }
      else
      {
        /* remove branch and merge adjacent inner nodes */
        pll_unode_t * prev = unode_prev(node);
        if (tree->vroot == node || tree->vroot == anode)
          tree->vroot = prev;

        /* find out position of to-be-removed node in the tree->nodes array,
         * and earmark it for deletion by setting respective entry to NULL */
        unsigned int anode_pos = clv2pos_map[anode->clv_index - tip_count];
        assert(anode_pos >= tip_count && anode_pos < node_count);
        tree->nodes[anode_pos] = NULL;
        tree->nodes[i] = prev;

        remove_branch(node);
        removed_count++;

        node = prev->next;
      }
    }
    while(node && node != start_node);
  }

  if (removed_count > 0)
  {
    /* compress tree->nodes array by excluding removed inner nodes */
    unsigned int idx = tip_count;
    unsigned int new_node_count = node_count - removed_count;
    for (unsigned int i = tip_count; i < node_count; ++i)
    {
      pll_unode_t * node  = tree->nodes[i];
      if (node)
        tree->nodes[idx++] = node;
    }
    assert(idx == new_node_count);

    /* update pll_utree_t metadata */
    tree->inner_count -= removed_count;
    tree->edge_count -= removed_count;
    tree->binary = 0;
    tree->nodes = (pll_unode_t **) realloc(tree->nodes,
                                           new_node_count*sizeof(pll_unode_t *));
  }

  free(clv2pos_map);

  return PLL_SUCCESS;
}


PLL_EXPORT int pllmod_utree_root_inplace(pll_utree_t * tree)
{
  if (!tree)
  {
    pllmod_set_error(PLL_ERROR_PARAM_INVALID, "Empty tree specified!");
    return PLL_FAILURE;
  }

  /* check if tree is already rooted */
  if (tree->vroot->next && tree->vroot->next->next == tree->vroot)
    return PLL_SUCCESS;

  pll_unode_t * root = tree->vroot;
  pll_unode_t * root_back = root->back;
  pll_unode_t * root_left = (pll_unode_t *) calloc(1, sizeof(pll_unode_t));
  pll_unode_t * root_right = (pll_unode_t *) calloc(1, sizeof(pll_unode_t));
  root_left->next = root_right;
  root_right->next = root_left;
  double root_brlen = root->length / 2.;
  unsigned int last_clv_index = 0;
  int last_scaler_index = 0;
  unsigned int last_node_index = 0;
  unsigned int last_pmatrix_index = 0;
  unsigned int node_count = tree->inner_count + tree->tip_count;

  for (unsigned int i = 0; i < node_count; ++i)
  {
    const pll_unode_t * node  = tree->nodes[i];
    last_clv_index = PLL_MAX(last_clv_index, node->clv_index);
    last_scaler_index = PLL_MAX(last_scaler_index, node->scaler_index);
    do
    {
      last_node_index = PLL_MAX(last_node_index, node->node_index);
      last_pmatrix_index = PLL_MAX(last_pmatrix_index, node->pmatrix_index);
      node = node->next;
    }
    while(node && node != tree->nodes[i]);
  }

  root_left->clv_index = root_right->clv_index = ++last_clv_index;
  root_left->scaler_index = root_right->scaler_index = ++last_scaler_index;
  root_left->node_index = ++last_node_index;
  root_right->node_index = ++last_node_index;
  root_right->pmatrix_index = ++last_pmatrix_index;

  pllmod_utree_connect_nodes(root, root_left, root_brlen);
  pllmod_utree_connect_nodes(root_right, root_back, root_brlen);

  tree->vroot = root_left;
  tree->inner_count++;
  tree->edge_count++;
  node_count++;

  tree->nodes = (pll_unode_t **) realloc(tree->nodes,
                                         node_count*sizeof(pll_unode_t *));
  tree->nodes[node_count-1] = root_left;

  return PLL_SUCCESS;
}

PLL_EXPORT int pllmod_utree_outgroup_root(pll_utree_t * tree,
                                          unsigned int * outgroup_tip_ids,
                                          unsigned int outgroup_size,
                                          int add_root_node)
{
  pll_unode_t ** split_to_node_map = NULL;
  pll_split_t * tree_splits = NULL;
  pll_unode_t * new_root = NULL;
  unsigned int tip_count;
  unsigned int split_count;

  if (!tree || !outgroup_tip_ids || !outgroup_size)
  {
    pllmod_set_error(PLL_ERROR_PARAM_INVALID,
                     "Empty tree and/or outgroup specified!");
    return PLL_FAILURE;
  }

  if (outgroup_size == 1)
  {
    // special case single-taxon outgroup: just find a tip by node_index
    for (unsigned int i = 0; i < tree->tip_count; ++i)
    {
      const pll_unode_t * node  = tree->nodes[i];
      if (node->node_index == outgroup_tip_ids[0])
      {
        new_root = node->back;
        break;
      }
    }
  }
  else
  {
    tip_count = tree->tip_count;
    split_count = tip_count - 3;

    split_to_node_map = (pll_unode_t **) calloc(split_count,
                                                sizeof(pll_unode_t *));

    if (!split_to_node_map)
    {
      pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                       "Cannot allocate memory for split->node map!");
      return PLL_FAILURE;
    }

    tree_splits = pllmod_utree_split_create(tree->vroot, tree->tip_count,
                                            split_to_node_map);

    if (!tree_splits)
    {
      assert(pll_errno);
      free(split_to_node_map);
      return PLL_FAILURE;
    }

    // create outgroup split
    pll_split_t outgroup_split = pllmod_utree_split_from_tips(outgroup_tip_ids,
                                                              outgroup_size,
                                                              tip_count);

    // check if this split is in the tree
    unsigned int split_len = bitv_length(tip_count);
    for (unsigned int i = 0; i < split_count; ++i)
    {
      if (!bitv_compare(tree_splits[i], outgroup_split, split_len))
      {
        new_root = split_to_node_map[i];
        break;
      }
    }

    pllmod_utree_split_destroy(tree_splits);
    free(split_to_node_map);
    free(outgroup_split);
  }

  // set tree->vroot to the outgroup split node
  if (new_root)
  {
    tree->vroot = new_root;
    if (add_root_node)
      return pllmod_utree_root_inplace(tree);
    else
      return PLL_SUCCESS;
  }
  else
  {
    pllmod_set_error(PLLMOD_TREE_ERROR_POLYPHYL_OUTGROUP,
                     "Outgroup is not monophyletic!");
    return PLL_FAILURE;
  }
}

int utree_insert_tips_random(pll_unode_t ** nodes, unsigned int taxa_count,
                             unsigned int start_tip, unsigned int random_seed)
{
  unsigned int i;
  unsigned int start_inner_count     = start_tip - 2;
  unsigned int start_branches        = 2 * start_tip - 3;
  unsigned int max_branches          = 2 * taxa_count - 3;
  unsigned int placed_branches_count = 0;
  unsigned int last_branch_id        = 0;

  pll_unode_t ** branches   = NULL;
  pll_random_state * rstate = NULL;

  branches = (pll_unode_t **) calloc(max_branches, sizeof(pll_unode_t *));

  if (!branches)
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for branches!");
    return PLL_FAILURE;
  }

  rstate =  pll_random_create(random_seed);

  if (!rstate)
  {
    free(branches);
    return PLL_FAILURE;
  }

  // check pmatrix indices on tip branches
  for (i = 0; i < taxa_count; ++i)
    last_branch_id = PLL_MAX(last_branch_id, nodes[i]->pmatrix_index);

  for (i = taxa_count; i < taxa_count + start_inner_count; ++i)
  {
    pll_unode_t * snode = nodes[i];
    do
    {
      if (snode->clv_index > snode->back->clv_index)
      {
        branches[placed_branches_count++] = snode;
        last_branch_id = PLL_MAX(last_branch_id, snode->pmatrix_index);
      }
      snode = snode->next;
    }
    while (snode != nodes[i]);
  }
  assert(placed_branches_count == start_branches);

  for (i = start_tip; i < taxa_count; ++i)
  {
    /* take tips iteratively */
    pll_unode_t * next_tip = nodes[i];
    pll_unode_t * next_inner = nodes[taxa_count + i - 2];

    /* select random branch from the tree */
    int rand_branch_id = pll_random_getint(rstate, placed_branches_count);
    pll_unode_t * next_branch = branches[rand_branch_id];

    /* connect tip to selected branch */
    pllmod_utree_connect_nodes(next_branch->back, next_inner,
                               PLLMOD_TREE_DEFAULT_BRANCH_LENGTH);
    pllmod_utree_connect_nodes(next_branch, next_inner->next,
                               PLLMOD_TREE_DEFAULT_BRANCH_LENGTH);
    pllmod_utree_connect_nodes(next_tip, next_inner->next->next,
                               PLLMOD_TREE_DEFAULT_BRANCH_LENGTH);

    if (pllmod_utree_is_tip (next_inner->back))
    {
      next_inner->next->pmatrix_index = next_inner->next->back->pmatrix_index =
          ++last_branch_id;
    }
    else
    {
      next_inner->pmatrix_index = next_inner->back->pmatrix_index =
          ++last_branch_id;
    }

    /* store branches */
    branches[placed_branches_count++] = next_inner;
    branches[placed_branches_count++] = next_inner->next->next;
  }
  assert(placed_branches_count == max_branches);

  /* clean */
  free (branches);
  pll_random_destroy(rstate);

  return PLL_SUCCESS;
}

/**
 * Extend a tree by inserting new taxa to randomly chosen branches
 */
PLL_EXPORT int pllmod_utree_extend_random(pll_utree_t * tree,
                                          unsigned int ext_taxa_count,
                                          const char * const* ext_names,
                                          unsigned int random_seed)
{
  unsigned int old_taxa_count  = tree->tip_count;
  unsigned int old_inner_count = tree->inner_count;
  unsigned int old_node_count  = old_taxa_count + old_inner_count;
  unsigned int new_taxa_count  = old_taxa_count + ext_taxa_count;
  unsigned int new_inner_count = old_inner_count + ext_taxa_count;
  unsigned int new_node_count = new_taxa_count + new_inner_count;

  unsigned int last_clv_id     = 0;
  unsigned int last_pmatrix_id = 0;
  unsigned int last_node_id    = 0;
  int next_scaler_id  = 0;

  unsigned int i;
  int retval;

  pll_unode_t ** old_nodes = tree->nodes;
  pll_unode_t ** new_nodes = (pll_unode_t **) calloc(new_node_count,
                                                     sizeof(pll_unode_t *));

  if (!new_nodes)
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for nodes!");
    return PLL_FAILURE;
  }

  // 1:1 mapping for old tips
  for (i = 0; i < old_taxa_count; ++i)
    new_nodes[i] = old_nodes[i];

  // copy old inner nodes and adjust clvs
  for (i = old_taxa_count; i < old_node_count; ++i)
  {
    unsigned int new_idx = i + ext_taxa_count;
    new_nodes[new_idx] = old_nodes[i];
    pll_unode_t * snode = new_nodes[new_idx];
    assert(snode->next);
    do
    {
      snode->clv_index += ext_taxa_count;
      snode->node_index += ext_taxa_count;
      last_clv_id = PLL_MAX(last_clv_id, snode->clv_index);
      last_node_id = PLL_MAX(last_node_id, snode->node_index);
      next_scaler_id = PLL_MAX(next_scaler_id, snode->scaler_index);
      last_pmatrix_id = PLL_MAX(last_pmatrix_id, snode->pmatrix_index);
      snode = snode->next;
    }
    while (snode != new_nodes[new_idx]);
  }

  // create new tip nodes
  for (i = old_taxa_count; i < new_taxa_count; ++i)
  {
    pll_unode_t * node = (pll_unode_t *) calloc(1, sizeof(pll_unode_t));
    node->clv_index = i;
    node->node_index = i;
    node->scaler_index = PLL_SCALE_BUFFER_NONE;
    node->pmatrix_index = ++last_pmatrix_id; // ????

    node->label = ext_names ? strdup(ext_names[i - old_taxa_count]) : NULL;

    new_nodes[i] = node;
  }

  // create new inner nodes
  for (i = old_node_count + ext_taxa_count; i < new_node_count; ++i)
  {
    pll_unode_t * node = pllmod_utree_create_node(++last_clv_id,
                                                  ++next_scaler_id,
                                                  NULL, NULL);

    node->node_index = ++last_node_id;
    node->next->node_index = ++last_node_id;
    node->next->next->node_index = ++last_node_id;

    new_nodes[i] = node;
  }

  retval = utree_insert_tips_random(new_nodes, new_taxa_count,
                                    old_taxa_count, random_seed);

  if (retval)
  {
    free(tree->nodes);
    tree->nodes = new_nodes;
    tree->tip_count = new_taxa_count;
    tree->inner_count = new_inner_count;
    tree->edge_count += 2 * ext_taxa_count;
    return PLL_SUCCESS;
  }
  else
  {
    free(new_nodes);
    return PLL_FAILURE;
  }
}

/**
 * Creates a random topology with default branch lengths
 */
PLL_EXPORT pll_utree_t * pllmod_utree_create_random(unsigned int taxa_count,
                                                    const char * const* names,
                                                    unsigned int random_seed)
{
  /*
   * The algorithm works as follows:
   *    1. Build a minimal 3-tip tree
   *    2. Select a branch at random
   *    3. Connect next tip to that branch
   *    4. Repeat 2 and 3 until no tips left
   */
  unsigned int i;
  unsigned int tip_node_count        = taxa_count;
  unsigned int inner_node_count      = taxa_count - 2;
  unsigned int node_count            = tip_node_count + inner_node_count;

  pll_unode_t ** nodes    = (pll_unode_t **) calloc(node_count,
                                                    sizeof(pll_unode_t *));

  pll_unode_t * tree_root;

  pll_utree_t * wrapped_tree;

  unsigned int node_id = 0;

  /* allocate tips */
  for (i=0; i<taxa_count; ++i)
  {
    nodes[i] = (pll_unode_t *)calloc(1, sizeof(pll_unode_t));
    nodes[i]->clv_index = i;
    nodes[i]->scaler_index = PLL_SCALE_BUFFER_NONE;
    nodes[i]->pmatrix_index = i;
    nodes[i]->node_index = node_id++;

    if (names)
    {
      nodes[i]->label = (char *) malloc( strlen(names[i]) + 1 );
      strcpy(nodes[i]->label, names[i]);
    }
    else
    {
      nodes[i]->label = NULL;
    }
  }

  /* allocate inner */
  for (i=taxa_count; i<node_count; ++i)
  {
    nodes[i] = pllmod_utree_create_node(i, (int)i, NULL, NULL);
    nodes[i]->scaler_index -= taxa_count;
    nodes[i]->next->scaler_index -= taxa_count;
    nodes[i]->next->next->scaler_index -= taxa_count;

    nodes[i]->node_index = node_id++;
    nodes[i]->next->node_index = node_id++;
    nodes[i]->next->next->node_index = node_id++;
  }
  assert(node_id == tip_node_count + inner_node_count * 3);

  /* set an inner node as return value */
  tree_root = nodes[taxa_count];

  /* build minimal tree with 3 tips and 1 inner node */
  pllmod_utree_connect_nodes(nodes[0], nodes[taxa_count],
                             PLLMOD_TREE_DEFAULT_BRANCH_LENGTH);
  pllmod_utree_connect_nodes(nodes[1], nodes[taxa_count]->next,
                             PLLMOD_TREE_DEFAULT_BRANCH_LENGTH);
  pllmod_utree_connect_nodes(nodes[2], nodes[taxa_count]->next->next,
                             PLLMOD_TREE_DEFAULT_BRANCH_LENGTH);

  /* insert remaining taxa_count-3 tips into the tree */
  utree_insert_tips_random(nodes, taxa_count, 3, random_seed);

  /* clean */
  free (nodes);

  wrapped_tree = pll_utree_wraptree(tree_root, tip_node_count);
  return (wrapped_tree);
}

/**
 * Creates a maximum parsimony topology using randomized stepwise-addition
 * algorithm. All branch lengths will be set to default.
 */
PLL_EXPORT
pll_utree_t * pllmod_utree_create_parsimony(unsigned int taxon_count,
                                            unsigned int seq_length,
                                            char * const * names,
                                            char * const * sequences,
                                            const unsigned int * site_weights,
                                            const pll_state_t * map,
                                            unsigned int states,
                                            unsigned int attributes,
                                            unsigned int random_seed,
                                            unsigned int * score)
{
  size_t i;
  pll_utree_t * tree = NULL;

  pll_partition_t * partition = pll_partition_create(taxon_count,
                                                     0,   /* number of CLVs */
                                                     states,
                                                     seq_length,
                                                     1,
                                                     1, /* pmatrix count */
                                                     1,  /* rate_cats */
                                                     0,  /* scale buffers */
                                                     attributes);

  if (!partition)
  {
    assert(pll_errno);
    return NULL;
  }

  /* set pattern weights and free the weights array */
  if (site_weights)
    pll_set_pattern_weights(partition, site_weights);

  /* find sequences in hash table and link them with the corresponding taxa */
  for (i = 0; i < taxon_count; ++i)
    pll_set_tip_states(partition, i, map, sequences[i]);

  tree = pllmod_utree_create_parsimony_multipart(taxon_count,
                                                 names,
                                                 1,
                                                 &partition,
                                                 random_seed,
                                                 score);

  /* destroy all structures allocated for the concrete PLL partition instance */
  pll_partition_destroy(partition);

  return tree;
}

/**
 * Creates a maximum parsimony topology using randomized stepwise-addition
 * algorithm. All branch lengths will be set to default.
 * This function can be used with partitioned alignments (e.g., combined DNA+AA data)
 */
PLL_EXPORT
pll_utree_t * pllmod_utree_create_parsimony_multipart(unsigned int taxon_count,
                                                      char * const * taxon_names,
                                                      unsigned int partition_count,
                                                      pll_partition_t * const * partitions,
                                                      unsigned int random_seed,
                                                      unsigned int * score)
{
  pll_utree_t * tree = NULL;
  unsigned int i;

  pll_parsimony_t ** parsimony =
      (pll_parsimony_t **) calloc(partition_count, sizeof(pll_parsimony_t *));

  if (!parsimony)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return NULL;
  }

  for (i = 0; i < partition_count; ++i)
  {
    assert(taxon_count == partitions[i]->tips);
    parsimony[i] = pll_fastparsimony_init(partitions[i]);
    if (!parsimony[i])
    {
      assert(pll_errno);
      goto cleanup;
    }
  }

  tree = pll_fastparsimony_stepwise(parsimony,
                                    taxon_names,
                                    score,
                                    partition_count,
                                    random_seed);

  if (tree)
  {
    /* update pmatrix/scaler/node indices */
    pll_utree_reset_template_indices(tree->nodes[tree->tip_count +
                                                 tree->inner_count - 1],
                                     tree->tip_count);

    /* set default branch lengths */
    pllmod_utree_set_length_recursive(tree,
                                      PLLMOD_TREE_DEFAULT_BRANCH_LENGTH,
                                      0);
  }
  else
    assert(pll_errno);

cleanup:
  /* destroy parsimony */
  for (i = 0; i < partition_count; ++i)
  {
    if (parsimony[i])
      pll_parsimony_destroy(parsimony[i]);
  }

  free(parsimony);

  return tree;
}

/* static functions */

static int utree_find_node_in_subtree(pll_unode_t * root,
                                                 pll_unode_t * node)
{
  if (root == node)
  {
    return PLL_SUCCESS;
  }

  if (root->next)
  {
    if (root->next == node || root->next->next == node)
    {
      return PLL_SUCCESS;
    }

    return utree_find_node_in_subtree(root->next->back, node)
        || utree_find_node_in_subtree(root->next->next->back, node);
  }

  return PLL_FAILURE;
}



/******************************************************************************/
/* Additional utilities */

static int utree_traverse_apply(pll_unode_t * node,
                                int (*cb_pre_trav)(pll_unode_t *, void *),
                                int (*cb_in_trav)(pll_unode_t *, void *),
                                int (*cb_post_trav)(pll_unode_t *, void *),
                                void *data)
{
  int retval = 1;
  pll_unode_t * child_tree = 0;

  if (cb_pre_trav && !cb_pre_trav(node,  data))
    return PLL_FAILURE;

  if (pllmod_utree_is_tip(node))
  {
    if (cb_in_trav)
      retval &= cb_in_trav(node, data);
    if (cb_post_trav)
      retval &= cb_post_trav(node, data);
    return retval;
  }

  child_tree = node->next;
  while(child_tree != node)
  {
    retval &= utree_traverse_apply(child_tree->back,
                                   cb_pre_trav, cb_in_trav, cb_post_trav, data);

    if (cb_in_trav &&
        child_tree->next != node &&
        !cb_in_trav(child_tree, data))
      return PLL_FAILURE;

    child_tree = child_tree->next;
  }

  if (cb_post_trav)
    retval &= cb_post_trav(node,  data);

  return retval;
}

PLL_EXPORT int pllmod_utree_traverse_apply(pll_unode_t * root,
                                           int (*cb_pre_trav)(pll_unode_t *,
                                                               void *),
                                           int (*cb_in_trav)(pll_unode_t *,
                                                               void *),
                                           int (*cb_post_trav)(pll_unode_t *,
                                                               void *),
                                           void *data)
{
  int retval = 1;

  assert(root);

  if (pllmod_utree_is_tip(root)) return PLL_FAILURE;

  retval &= utree_traverse_apply(root->back,
                                 cb_pre_trav, cb_in_trav, cb_post_trav,
                                 data);
  retval &= utree_traverse_apply(root,
                                 cb_pre_trav, cb_in_trav, cb_post_trav,
                                 data);

  return retval;
}

PLL_EXPORT int pllmod_utree_is_tip(const pll_unode_t * node)
{
  return (node->next == NULL);
}

PLL_EXPORT int pllmod_rtree_is_tip(const pll_rnode_t * node)
{
  return (node->left == NULL && node->right == NULL);
}

PLL_EXPORT void pllmod_utree_set_length(pll_unode_t * edge,
                                            double length)
{
  edge->length = edge->back->length = length;
}

PLL_EXPORT void pllmod_utree_set_length_recursive(pll_utree_t * tree,
                                                  double length,
                                                  int missing_only)
{
  /* set branch lengths */
  unsigned int i;
  unsigned int tip_count = tree->tip_count;
  unsigned int inner_count = tree->inner_count;
  for (i = 0; i < tip_count + inner_count; ++i)
  {
    pll_unode_t * node = tree->nodes[i];
    if (!node->length || !missing_only)
      pllmod_utree_set_length(node, length);
    if (node->next)
    {
      if (!node->next->length || !missing_only)
        pllmod_utree_set_length(node->next, length);
      if (!node->next->next->length || !missing_only)
        pllmod_utree_set_length(node->next->next, length);
    }
  }
}

PLL_EXPORT void pllmod_utree_scale_branches(pll_utree_t * tree,
                                            double branch_length_scaler)
{
  /* scale branch lengths */
  unsigned int i;
  unsigned int tip_count = tree->tip_count;
  unsigned int inner_count = tree->inner_count;
  pll_unode_t ** nodes = tree->nodes;
  for (i=0; i<tip_count; ++i)
  {
    nodes[i]->length *= branch_length_scaler;
  }
  for (i=tip_count; i<tip_count+inner_count; ++i)
  {
    nodes[i]->length *= branch_length_scaler;
    nodes[i]->next->length *= branch_length_scaler;
    nodes[i]->next->next->length *= branch_length_scaler;
  }
}

PLL_EXPORT void pllmod_utree_scale_branches_all(pll_unode_t * root,
                                                double branch_length_scaler)
{
  double root_length;

  /* scale all branches in a tree */
  pllmod_utree_scale_subtree_branches(root, branch_length_scaler);
  root_length = root->length;
  pllmod_utree_scale_subtree_branches(root->back, branch_length_scaler);

  /* undo duplicated scaling */
  root->length = root->back->length = root_length;
}

PLL_EXPORT void pllmod_utree_scale_subtree_branches(pll_unode_t * root,
                                                    double branch_length_scaler)
{
 /* scale all branches in a subtree rooted at node */
  root->length *= branch_length_scaler;
  root->back->length *= branch_length_scaler;

  if (root->next)
  {
    pllmod_utree_scale_subtree_branches(root->next->back, branch_length_scaler);
    pllmod_utree_scale_subtree_branches(root->next->next->back, branch_length_scaler);
  }
}


/**
 * compute the likelihood on a utree structure
 * if update_pmatrices or update_partials are set, p-matrices and CLVs are
 * updated before computing the likelihood.
 */
PLL_EXPORT double pllmod_utree_compute_lk(pll_partition_t * partition,
                                       pll_unode_t * tree,
                                       const unsigned int * params_indices,
                                       int update_pmatrices,
                                       int update_partials)
{
  struct cb_params parameters;
  assert (tree);
  assert (tree->pmatrix_index == tree->back->pmatrix_index);

  parameters.partition      = partition;
  parameters.params_indices = params_indices;

  /* update pmatrices */
  if (update_pmatrices || update_partials)
  {
    parameters.update_pmatrices = update_pmatrices;
    parameters.update_partials  = update_partials;

    pllmod_utree_traverse_apply(tree,
                                0,
                                0,
                                cb_update_matrices_partials,
                                (void *) &parameters);
  }

  double logl = pll_compute_edge_loglikelihood(partition,
                                              tree->clv_index,
                                              tree->scaler_index,
                                              tree->back->clv_index,
                                              tree->back->scaler_index,
                                              tree->pmatrix_index,
                                              params_indices,
                                              NULL);
  return logl;
}

struct clv_set_data
{
  int * set_indices;
  unsigned int max_index;
  unsigned int tip_count;
};

static int cb_set_clv_minimal(pll_unode_t * node, void * data)
{
  unsigned int i, next_index;
  int index_found;
  struct clv_set_data * clv_data = (struct clv_set_data *)data;
  int * v = 0;

  if (!pllmod_utree_is_tip(node))
  {
    /* find next free position */
    v = clv_data->set_indices;
    next_index  = 0;
    index_found = 0;
    for (i=0; i<clv_data->max_index; ++i)
    {
      if (!v[i])
      {
        index_found = 1;
        next_index = i;
        v[i] = 1;
        break;
      }
    }
    assert(index_found);

    /* set clv index */
    node->clv_index =
      node->next->clv_index =
      node->next->next->clv_index =
       next_index + clv_data->tip_count;
    /* set scaler index */
    node->scaler_index =
       node->next->scaler_index =
       node->next->next->scaler_index =
        (int)(next_index + clv_data->tip_count);

    /* free indices from children */
    if (!pllmod_utree_is_tip(node->next->back))
    {
      v[node->next->back->clv_index - clv_data->tip_count] = 0;
    }
    if (!pllmod_utree_is_tip(node->next->next->back))
    {
      v[node->next->next->back->clv_index - clv_data->tip_count] = 0;
    }
  }

  /* continue */
  return 1;
}

PLL_EXPORT int pllmod_utree_set_clv_minimal(pll_unode_t * root,
                                            unsigned int tip_count)
{
  unsigned int clv_count = (unsigned int) ceil(log2(tip_count)) + 2;
  int * set_indices = (int *) calloc((size_t)clv_count, sizeof(int));
  struct clv_set_data data;
  data.set_indices = set_indices;
  data.max_index   = clv_count;
  data.tip_count   = tip_count;
  pllmod_utree_traverse_apply(root, 0, 0, cb_set_clv_minimal, (void *) &data);
  free(set_indices);

  return PLL_SUCCESS;
}

static int rtree_traverse_apply(pll_rnode_t * node,
                                int (*cb_pre_trav)(pll_rnode_t *, void *),
                                int (*cb_in_trav)(pll_rnode_t *, void *),
                                int (*cb_post_trav)(pll_rnode_t *, void *),
                                void *data)
{
  int retval = 1;

  if (cb_pre_trav && !cb_pre_trav(node,  data))
    return PLL_FAILURE;

  if (node->left)
  {
    retval &= rtree_traverse_apply(node->left,
                                   cb_pre_trav,
                                   cb_in_trav,
                                   cb_post_trav,
                                   data);
  }

  if (cb_in_trav && !cb_in_trav(node,  data))
    return PLL_FAILURE;

  if (node->right) {
    retval &= rtree_traverse_apply(node->right,
                                   cb_pre_trav,
                                   cb_in_trav,
                                   cb_post_trav,
                                   data);
  }

  if (cb_post_trav)
    retval &= cb_post_trav(node,  data);

  return retval;
}

PLL_EXPORT int pllmod_rtree_traverse_apply(pll_rnode_t * root,
                                           int (*cb_pre_trav)(pll_rnode_t *,
                                                               void *),
                                           int (*cb_in_trav)(pll_rnode_t *,
                                                               void *),
                                           int (*cb_post_trav)(pll_rnode_t *,
                                                               void *),
                                           void *data)
{
  int retval = 1;

  if (!root->left || !root->right) return PLL_FAILURE;

  retval &= rtree_traverse_apply(root,
                                 cb_pre_trav,
                                 cb_in_trav,
                                 cb_post_trav,
                                 data);

  return retval;
}

/* auxiliary structure for the callback function below */
struct serial_tree_s {
  pll_unode_t * serialized_tree;
  unsigned int node_count;
  unsigned int max_nodes;
};

/* callback function to fill the serialized tree */
static int cb_serialize(pll_unode_t * tree,
                        void * data)
{
  struct serial_tree_s * list = (struct serial_tree_s *) data;
  pll_unode_t * serialized_tree = list->serialized_tree;
  unsigned int cur_pos = list->node_count;

  assert(cur_pos < list->max_nodes);

  memcpy(&(serialized_tree[cur_pos]), tree, sizeof(pll_unode_t));
  serialized_tree[cur_pos].data  = 0;
  serialized_tree[cur_pos].label = 0;
  if (!pllmod_utree_is_tip(tree))
  {
    /* set to arbitrary non-junk value */
    serialized_tree[cur_pos].next = (pll_unode_t *) 1;
  }

  ++list->node_count;
  return 1;
}

//TODO: serialize/expand using a compressed format instead of pll_unode_t
PLL_EXPORT pll_unode_t * pllmod_utree_serialize(pll_unode_t * tree,
                                                unsigned int tip_count)
{
  unsigned int node_count;
  pll_unode_t * serialized_tree;
  struct serial_tree_s data;

  node_count = 2*tip_count - 2;

  /* allocate the serialized structure */
  serialized_tree = (pll_unode_t *) malloc(node_count * sizeof (pll_unode_t));

  /* fill data for callback function */
  data.serialized_tree = serialized_tree;
  data.node_count      = 0;
  data.max_nodes       = node_count;

  /* if tree is a tip, move to its back position */
  if (pllmod_utree_is_tip(tree)) tree = tree->back;

  /* apply callback function to serialize */
  pllmod_utree_traverse_apply(tree,
                              NULL,
                              NULL,
                              cb_serialize,
                              &data);

  if (data.node_count != data.max_nodes)
  {
    /* if the number of serialized nodes is not correct, return error */
    pllmod_set_error(PLLMOD_TREE_ERROR_INVALID_TREE,
                     "tree structure ot tip_count are invalid");
    free(serialized_tree);
    serialized_tree = NULL;
  }

  return serialized_tree;
}

PLL_EXPORT pll_utree_t * pllmod_utree_expand(pll_unode_t * serialized_tree,
                                             unsigned int tip_count)
{
  unsigned int i, node_count, next_node_index;
  pll_unode_t ** tree_stack;
  pll_unode_t * tree;
  unsigned int tree_stack_top;

  pllmod_reset_error();

  node_count  = 2*tip_count - 2;

  /* allocate stack for at most 'n_tips' nodes */
  tree_stack = (pll_unode_t **) malloc(tip_count * sizeof (pll_unode_t *));
  tree_stack_top = 0;

  next_node_index = tip_count;

  /* read nodes */
  for (i=0; i<node_count; ++i)
  {
    pll_unode_t * t = 0;                   /* new node */
    pll_unode_t t_s = serialized_tree[i];  /* serialized node */
    if (t_s.next)
    {
      /* build inner node and connect */
      pll_unode_t *t_cr, *t_r, *t_cl, *t_l;
      t = pllmod_utree_create_node(t_s.clv_index,
                                   t_s.scaler_index,
                                   0,  /* label */
                                   0); /* data */
      t_l = t->next;
      t_r = t->next->next;

      t->node_index = next_node_index++;
      t_r->node_index = next_node_index++;
      t_l->node_index = next_node_index++;

      /* pop and connect */
      t_cr = tree_stack[--tree_stack_top];
      t_r->back = t_cr; t_cr->back = t_r;
      t_cl = tree_stack[--tree_stack_top];
      t_l->back = t_cl; t_cl->back = t_l;

      /* set branch attributes */
      t->length = t_s.length;
      t->pmatrix_index = t_s.pmatrix_index;
      t_r->pmatrix_index = t_cr->pmatrix_index;
      t_r->length = t_cr->length;
      t_l->pmatrix_index = t_cl->pmatrix_index;
      t_l->length = t_cl->length;
    }
    else
    {
      t = (pll_unode_t *)calloc(1, sizeof(pll_unode_t));
      memcpy(t, &t_s, sizeof(pll_unode_t));
      assert(t->node_index < tip_count);
    }

    /* push */
    tree_stack[tree_stack_top++] = t;
  }

  /* root vertices must be in the stack */
  assert (tree_stack_top == 2);
  assert (next_node_index == (4*tip_count - 6));

  tree = tree_stack[--tree_stack_top];
  tree->back = tree_stack[--tree_stack_top];
  tree->back->back = tree;

  if(tree->pmatrix_index != tree->back->pmatrix_index)
  {
    /* if pmatrix indices differ, connecting branch must be a tip */
    if(tree->back->next)
    {
      pllmod_set_error(PLLMOD_TREE_ERROR_INVALID_TREE,
                       "pmatrix indices do not match in serialized tree");
    }
    tree->pmatrix_index = tree->back->pmatrix_index;
  }

  if(tree->length != tree->back->length)
  {
    pllmod_set_error(PLLMOD_TREE_ERROR_INVALID_TREE,
                     "branch lengths do not matchin serialized tree");
  }

  if (pll_errno)
  {
    pll_utree_graph_destroy(tree, NULL);
    tree = 0;
  }

  free(tree_stack);

  return pll_utree_wraptree(tree, tip_count);
}

PLL_EXPORT int pllmod_utree_draw_support(pll_utree_t * ref_tree,
                                         const double * support,
                                         pll_unode_t ** node_map,
                                         char * (*cb_serialize)(double) )
{
  if (!ref_tree || !support)
  {
    pllmod_set_error(PLL_ERROR_PARAM_INVALID, "Parameter is NULL!\n");
    return PLL_FAILURE;
  }

  unsigned int split_count = ref_tree->edge_count - ref_tree->tip_count;
  for (size_t i = 0; i < split_count; ++i)
  {
    pll_unode_t * node = node_map ? node_map[i] :
                                    ref_tree->nodes[ref_tree->tip_count + i];

    /* this has to be an inner node! */
    assert(node->next);

    if (node->label)
      free(node->label);

    node->label = node->next->label = node->next->next->label =
        cb_serialize ? cb_serialize(support[i]) : default_support_fmt(support[i]);
  }

  return PLL_SUCCESS;
}


/******************************************************************************/
/* Static functions */

static int rtree_rollback_tbr(pll_tree_rollback_t * rollback_info)
{
  PLLMOD_UNUSED(rollback_info);
  assert(UNIMPLEMENTED);
  return PLL_FAILURE;
}

static int rtree_rollback_spr(pll_tree_rollback_t * rollback_info)
{
  //TODO: Add preconditions

  pll_rnode_t * p = (pll_rnode_t *) rollback_info->SPR.prune_edge;
  pll_rnode_t * r = (pll_rnode_t *) rollback_info->SPR.regraft_edge;

  /* undo move */
  if (!pllmod_rtree_spr(p, r, 0, 0))
    return PLL_FAILURE;

  //TODO: set branch lengths
  //
  return PLL_SUCCESS;
}

static int rtree_rollback_nni(pll_tree_rollback_t * rollback_info)
{
  PLLMOD_UNUSED(rollback_info);
  assert(UNIMPLEMENTED);
  return PLL_FAILURE;
}


static int utree_rollback_tbr(pll_tree_rollback_t * rollback_info)
{
  assert(!rollback_info->rooted);
  assert(rollback_info->rearrange_type == PLLMOD_TREE_REARRANGE_TBR);

  pll_unode_t * p = (pll_unode_t *) rollback_info->TBR.bisect_edge;
  pll_unode_t * q = p->next->back;
  pll_unode_t * r = p->back->next->back;
  double reconn_length = rollback_info->TBR.reconn_edge.length;

  /* undo move */
    if (!pllmod_utree_tbr(p, &(rollback_info->TBR.reconn_edge), 0))
      return PLL_FAILURE;

  /* reset branches */
  pllmod_utree_set_length(p, reconn_length);
  pllmod_utree_set_length(q, rollback_info->TBR.bisect_left_bl);
  pllmod_utree_set_length(r, rollback_info->TBR.bisect_right_bl);
  pllmod_utree_set_length(p->next, rollback_info->TBR.reconn_parent_left_bl);
  pllmod_utree_set_length(p->next->next, rollback_info->TBR.reconn_parent_right_bl);
  pllmod_utree_set_length(p->back->next, rollback_info->TBR.reconn_child_left_bl);
  pllmod_utree_set_length(p->back->next->next, rollback_info->TBR.reconn_child_right_bl);

  return PLL_SUCCESS;
}

static int utree_rollback_spr(pll_tree_rollback_t * rollback_info)
{
  assert(!rollback_info->rooted);
  assert(rollback_info->rearrange_type == PLLMOD_TREE_REARRANGE_SPR);

  pll_unode_t * p = (pll_unode_t *) rollback_info->SPR.prune_edge;
  pll_unode_t * r = (pll_unode_t *) rollback_info->SPR.regraft_edge;
  pll_unode_t * z1 =  p->next->back;
  pll_unode_t * z2 =  r->back;

  /* undo move */
  if (!pllmod_utree_spr(p, r, 0))
    return PLL_FAILURE;

  /* reset branches */
  pllmod_utree_set_length(z1, rollback_info->SPR.regraft_bl);
  pllmod_utree_set_length(p, rollback_info->SPR.prune_bl);
  pllmod_utree_set_length(r, rollback_info->SPR.prune_left_bl);
  pllmod_utree_set_length(z2, rollback_info->SPR.prune_right_bl);

  return PLL_SUCCESS;
}

static int utree_rollback_nni(pll_tree_rollback_t * rollback_info)
{
  assert(!rollback_info->rooted);
  assert(rollback_info->rearrange_type == PLLMOD_TREE_REARRANGE_NNI);

  pll_unode_t * p = rollback_info->NNI.edge;
  pll_unode_t * q = p->back;

  /* undo move */
  if (!pllmod_utree_nni(p, rollback_info->NNI.type, 0))
      return PLL_FAILURE;

  /* reset branches */

  pllmod_utree_set_length(p, rollback_info->NNI.edge_bl);
  pllmod_utree_set_length(p->next, rollback_info->NNI.left_left_bl);
  pllmod_utree_set_length(p->next->next, rollback_info->NNI.left_right_bl);
  pllmod_utree_set_length(q->next, rollback_info->NNI.right_left_bl);
  pllmod_utree_set_length(q->next->next, rollback_info->NNI.right_right_bl);

  //assert(UNIMPLEMENTED);
  return PLL_SUCCESS;
}

/**
 * callback function for updating p-matrices and partials
 */
static int cb_update_matrices_partials(pll_unode_t * node, void *data)
{
  struct cb_params * st_data = (struct cb_params *) data;

  if (st_data->update_pmatrices)
  {
    unsigned int matrix_index = node->pmatrix_index;
    double branch_length = node->length;

    /* check integrity */
    assert(fabs(node->length - node->back->length) < 1e-8);
    assert(node->pmatrix_index == node->back->pmatrix_index);

    pll_update_prob_matrices (st_data->partition,
                              st_data->params_indices,
                              &matrix_index,
                              &branch_length,
                              1);
  }

  if (st_data->update_partials && !pllmod_utree_is_tip(node))
  {
    /* check integrity */
    assert(node->next->pmatrix_index == node->next->back->pmatrix_index);
    assert(node->next->next->pmatrix_index == node->next->next->back->pmatrix_index);

    pll_operation_t op;
    op.child1_clv_index    = node->next->back->clv_index;
    op.child1_scaler_index = node->next->back->scaler_index;
    op.child1_matrix_index = node->next->pmatrix_index;
    op.child2_clv_index    = node->next->next->back->clv_index;
    op.child2_scaler_index = node->next->next->back->scaler_index;
    op.child2_matrix_index = node->next->next->pmatrix_index;
    op.parent_clv_index    = node->clv_index;
    op.parent_scaler_index = node->scaler_index;

    pll_update_partials(st_data->partition, &op, 1);
  }

  return PLL_SUCCESS;
}

static void shuffle_tree_nodes(const pll_utree_t * tree, unsigned int seed)
{
  unsigned int node_count = tree->tip_count + tree->inner_count;
  pll_random_state * rstate =  pll_random_create(seed);
  pll_unode_t ** subnodes =  (pll_unode_t **) calloc(tree->tip_count,
                                                     sizeof(pll_unode_t *));

  for (unsigned int i = tree->tip_count; i < node_count; ++i)
  {
    pll_unode_t * node = tree->nodes[i];
    unsigned int degree = 0;
    do
    {
      subnodes[degree] = node;
      degree++;
      node = node->next;
    }
    while (node != tree->nodes[i]);

    // Fisher–Yates shuffle
    for (unsigned int j = degree-1; j > 0; --j)
    {
      unsigned int r = pll_random_getint(rstate, j+1);
      PLL_SWAP(subnodes[j], subnodes[r]);
    }

    // re-connect pll_unodes in the new, shuffled order
    tree->nodes[i] = node = subnodes[0];
    for (unsigned int j = 1; j < degree; ++j)
    {
      node->next = subnodes[j];
      node = node->next;
    }

    // close roundabout
    node->next = tree->nodes[i];
  }

  pll_random_destroy(rstate);
  free(subnodes);
}


static void split_multi_node(pll_utree_t * tree, pll_unode_t * first,
                             pll_unode_t * last, unsigned int degree)
{
  assert(last->next == first);
  if (degree > 3)
  {
    assert(first->next && first->next->next && first->next->next->next != first);

    // got a multifurcating node, split it in two
    unsigned int new_pmatrix_id = tree->edge_count;
    unsigned int new_node_id = tree->edge_count * 2;
    unsigned int new_clv_id = tree->tip_count + tree->inner_count;
    unsigned int new_scaler_id = tree->inner_count;

    pll_unode_t * second = first->next;

    pll_unode_t * old_link = (pll_unode_t *) calloc(1, sizeof(pll_unode_t));
    pll_unode_t * new_link = (pll_unode_t *) calloc(1, sizeof(pll_unode_t));

    old_link->data = second->data;
    new_link->data = second->next->data;

    //close 'new' roundabout
    new_link->next = second->next;
    last->next = new_link;

    //close 'old' roundabout
    old_link->next = first;
    second->next = old_link;

    old_link->clv_index = second->clv_index;
    old_link->scaler_index = second->scaler_index;
    old_link->pmatrix_index = new_pmatrix_id;
    old_link->node_index = new_node_id;

    assert(new_link->next && new_link->next->next);

    new_link->pmatrix_index = new_pmatrix_id;
    new_link->node_index = new_node_id+1;

    new_link->clv_index = new_link->next->clv_index =
        new_link->next->next->clv_index = new_clv_id;
    new_link->scaler_index = new_link->next->scaler_index =
        new_link->next->next->scaler_index = (int) new_scaler_id;

    //set backpointers old<->new
    pllmod_utree_connect_nodes(old_link, new_link, PLLMOD_TREE_DEFAULT_BRANCH_LENGTH);

    tree->nodes[tree->inner_count + tree->tip_count] = new_link;

    tree->edge_count++;
    tree->inner_count++;

    // split new node if needed
    split_multi_node(tree, new_link, last, degree-1);
  }
}

static char * default_support_fmt(double support)
{
  char *sup_str;
  int size_alloced = asprintf(&sup_str, "%lf", support);

  return size_alloced >= 0 ? sup_str : NULL;
}
