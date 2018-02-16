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

/**
 * Creates a random topology with default branch lengths
 */
PLL_EXPORT pll_utree_t * pllmod_utree_create_random(unsigned int taxa_count,
                                                    const char * const* names)
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
  unsigned int max_branches          = 2 * tip_node_count - 3;
  unsigned int placed_branches_count = 0;

  pll_unode_t ** nodes    = (pll_unode_t **) calloc(node_count,
                                                    sizeof(pll_unode_t *));
  pll_unode_t ** branches = (pll_unode_t **) calloc(max_branches,
                                                    sizeof(pll_unode_t *));

  pll_unode_t * next_tip;
  pll_unode_t * next_inner;
  pll_unode_t * next_branch;
  pll_unode_t * tree_root;

  pll_utree_t * wrapped_tree;

  unsigned int next_branch_id = taxa_count;
  unsigned int rand_branch_id;
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
  branches[placed_branches_count++] = nodes[taxa_count];
  pllmod_utree_connect_nodes(nodes[1], nodes[taxa_count]->next,
                             PLLMOD_TREE_DEFAULT_BRANCH_LENGTH);
  branches[placed_branches_count++] = nodes[taxa_count]->next;
  pllmod_utree_connect_nodes(nodes[2], nodes[taxa_count]->next->next,
                             PLLMOD_TREE_DEFAULT_BRANCH_LENGTH);
  branches[placed_branches_count++] = nodes[taxa_count]->next->next;

  for (i=3; i<taxa_count; ++i)
  {
    /* take tips iteratively */
    next_tip = nodes[i];
    next_inner = nodes[tip_node_count + i - 2];

    /* select random branch from the tree */
    rand_branch_id = (unsigned int) rand() % placed_branches_count;
    next_branch = branches[rand_branch_id];

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
          next_branch_id++;
    }
    else
    {
      next_inner->pmatrix_index = next_inner->back->pmatrix_index =
          next_branch_id++;
    }

    /* store branches */
    branches[placed_branches_count++] = next_inner;
    branches[placed_branches_count++] = next_inner->next->next;
  }
  assert(placed_branches_count == max_branches);

  /* clean */
  free (nodes);
  free (branches);

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

PLL_EXPORT int pllmod_utree_is_tip(pll_unode_t * node)
{
  return (node->next == NULL);
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
  pll_unode_t ** nodes = tree->nodes;
  for (i=0; i<tip_count; ++i)
  {
    if (!nodes[i]->length || !missing_only)
      pllmod_utree_set_length(nodes[i], length);
  }
  for (i=tip_count; i<tip_count+inner_count; ++i)
  {
    if (!nodes[i]->length || !missing_only)
      pllmod_utree_set_length(nodes[i], length);
    if (!nodes[i]->length || !missing_only)
      pllmod_utree_set_length(nodes[i], length);
    if (!nodes[i]->length || !missing_only)
      pllmod_utree_set_length(nodes[i], length);
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

    if (cb_in_trav && !cb_in_trav(node,  data))
      return PLL_FAILURE;

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
  int node_count;
  int max_nodes;
};

/* callback function to fill the serialized tree */
static int cb_serialize(pll_unode_t * tree,
                        void * data)
{
  struct serial_tree_s * list = (struct serial_tree_s *) data;
  pll_unode_t * serialized_tree = list->serialized_tree;
  int cur_pos = list->node_count;

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
