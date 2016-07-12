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
#include "pll_tree.h"
#include "../pllmod_common.h"

#define UNIMPLEMENTED 0

static int rtree_rollback_TBR(pll_tree_rollback_t * rollback_info);
static int rtree_rollback_SPR(pll_tree_rollback_t * rollback_info);
static int rtree_rollback_NNI(pll_tree_rollback_t * rollback_info);
static int utree_rollback_TBR(pll_tree_rollback_t * rollback_info);
static int utree_rollback_SPR(pll_tree_rollback_t * rollback_info);
static int utree_rollback_NNI(pll_tree_rollback_t * rollback_info);
static int utree_find_node_in_subtree(pll_utree_t * root,
                                                 pll_utree_t * node);
static int cb_update_matrices_partials(pll_utree_t * node, void *data);

struct cb_params
{
  unsigned int * params_indices;
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
 *
 * @returns true, if the move was applied correctly
 */
PLL_EXPORT int pll_utree_TBR(pll_utree_t * b_edge,
                             pll_tree_edge_t * r_edge,
                             pll_tree_rollback_t * rollback_info)
{
  pll_utree_t *parent, *child;

  /* validate if the move can be applied */

  /* 1. bisection point must not be a leaf branch */
  if (!(b_edge->next && b_edge->back->next))
  {
    pllmod_set_error(PLL_ERROR_TBR_LEAF_BISECTION, "attempting to bisect at a leaf node");
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
    pllmod_set_error(PLL_ERROR_TBR_OVERLAPPED_NODES, "TBR nodes are overlapped");
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
    pllmod_set_error(PLL_ERROR_TBR_SAME_SUBTREE, "TBR reconnection in same subtree");
    return PLL_FAILURE;
  }

  /* save rollback information */
  if (rollback_info)
  {
    rollback_info->rearrange_type     = PLL_REARRANGE_TBR;
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
  pll_utree_bisect(b_edge, &parent, &child);

  /* reconnect at r_edge */
  pll_utree_reconnect(r_edge,
                      b_edge);

  return PLL_SUCCESS;
}



/**
 * Performs one SPR move
 * The CLV, scaler and pmatrix indices are updated.
 *
 * @param[in] p_edge Edge to be pruned
 * @param[in] r_edge Edge to be regrafted
 *
 * @returns true, if the move was applied correctly
 */
PLL_EXPORT int pll_utree_SPR(pll_utree_t * p_edge,
                             pll_utree_t * r_edge,
                             pll_tree_rollback_t * rollback_info)
{
  int retval;
  pll_utree_t * pruned_tree;

  if (pll_utree_is_tip(p_edge))
  {
    /* invalid move */
    pllmod_set_error(PLL_ERROR_SPR_INVALID_NODE, "Attempting to prune a leaf branch");
    return PLL_FAILURE;
  }

  /* save rollback information */
  if (rollback_info)
  {
    rollback_info->rearrange_type     = PLL_REARRANGE_SPR;
    rollback_info->rooted             = 0;
    rollback_info->SPR.prune_edge     = (void *) p_edge;
    rollback_info->SPR.regraft_edge   = (void *) p_edge->next->back;
    rollback_info->SPR.prune_bl       = p_edge->length;
    rollback_info->SPR.prune_left_bl  = p_edge->next->length;
    rollback_info->SPR.prune_right_bl = p_edge->next->next->length;
    rollback_info->SPR.regraft_bl     = r_edge->length;
  }

  /* prune p_edge */
  pruned_tree = pll_utree_prune(p_edge);
  if (!pruned_tree)
  {
    /* check that errno was set correctly */
    assert(pll_errno & PLL_ERROR_SPR_MASK);
    return PLL_FAILURE;
  }

  /* regraft p_edge on r_edge*/
  retval = pll_utree_regraft(p_edge, r_edge);
  assert(retval == PLL_SUCCESS || pll_errno & PLL_ERROR_SPR_MASK);

  return retval;
}

/**
 * Performs one NNI move
 * The CLV, scaler and pmatrix indices are updated.
 *
 * @param[in] edge NNI interchange edge
 * @param[in] type move type: PLL_NNI_LEFT, PLL_NNI_RIGHT
 * @param[out] rollback_info Rollback information
 *
 * @returns true, if the move was applied correctly
 */
PLL_EXPORT int pll_utree_NNI(pll_utree_t * edge,
                             int type,
                             pll_tree_rollback_t * rollback_info)
{
  pll_utree_t *base, *neighbor;

  /* validate preconditions */
  assert(edge && edge->back);

  if (!(type == PLL_NNI_LEFT || type == PLL_NNI_RIGHT))
  {
    /* invalid move */
    pllmod_set_error(PLL_ERROR_NNI_INVALID_MOVE, "Invalid NNI move type");
    return PLL_FAILURE;
  }
  if (pll_utree_is_tip(edge) || pll_utree_is_tip(edge->back))
  {
    /* invalid move */
    pllmod_set_error(PLL_ERROR_INTERCHANGE_LEAF, "Attempting to apply NNI on a leaf branch");
    return PLL_FAILURE;
  }

  /* save rollback information */
  if (rollback_info)
  {
    rollback_info->rearrange_type     = PLL_REARRANGE_NNI;
    rollback_info->rooted             = 0;
    rollback_info->NNI.edge           = (void *) edge;
    rollback_info->NNI.type           = type;
    rollback_info->NNI.left_left_bl   = edge->next->length;
    rollback_info->NNI.left_right_bl  = edge->next->next->length;
    rollback_info->NNI.right_left_bl  = edge->back->next->length;
    rollback_info->NNI.right_right_bl = edge->back->next->next->length;
    rollback_info->NNI.edge_bl        = edge->length;
  }

  /* select node to interchange */
  base = edge->next;
  if (type == PLL_NNI_LEFT)
  {
    neighbor = edge->back->next;
  }
  else if (type == PLL_NNI_RIGHT)
  {
    neighbor = edge->back->next->next;
  }
  else
    assert(0);

  /* apply move */
  if (!pll_utree_interchange(base, neighbor))
    return PLL_FAILURE;

  return PLL_SUCCESS;
}

PLL_EXPORT int pll_tree_rollback(pll_tree_rollback_t * rollback_info)
{
  int retval = PLL_FAILURE;
  switch (rollback_info->rearrange_type)
  {
    case PLL_REARRANGE_TBR:
      {
        if (rollback_info->rooted)
          retval = rtree_rollback_TBR (rollback_info);
        else
          retval = utree_rollback_TBR (rollback_info);
      }
      break;
    case PLL_REARRANGE_SPR:
      {
        if (rollback_info->rooted)
          retval = rtree_rollback_SPR (rollback_info);
        else
          retval = utree_rollback_SPR (rollback_info);
      }
      break;
    case PLL_REARRANGE_NNI:
      {
        if (rollback_info->rooted)
          retval = rtree_rollback_NNI (rollback_info);
        else
          retval = utree_rollback_NNI (rollback_info);
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
PLL_EXPORT pll_utree_t * pll_utree_create_random(unsigned int n_taxa,
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
  unsigned int n_tip_nodes       = n_taxa;
  unsigned int n_inner_nodes     = n_taxa - 2;
  unsigned int n_nodes           = n_tip_nodes + n_inner_nodes;
  unsigned int max_branches      = 2 * n_tip_nodes - 3;
  unsigned int n_placed_branches = 0;

  pll_utree_t ** nodes    = (pll_utree_t **) calloc(n_nodes, sizeof(pll_utree_t *));
  pll_utree_t ** branches = (pll_utree_t **) calloc(max_branches, sizeof(pll_utree_t *));

  pll_utree_t * next_tip;
  pll_utree_t * next_inner;
  pll_utree_t * next_branch;
  pll_utree_t * new_tree;

  unsigned int next_branch_id = n_taxa;
  unsigned int rand_branch_id;

  /* allocate tips */
  for (i=0; i<n_taxa; ++i)
  {
    nodes[i] = (pll_utree_t *)calloc(1, sizeof(pll_utree_t));
    nodes[i]->clv_index = i;
    nodes[i]->scaler_index = PLL_SCALE_BUFFER_NONE;
    nodes[i]->pmatrix_index = i;

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
  for (i=n_taxa; i<n_nodes; ++i)
  {
    nodes[i] = pll_utree_create_node(i, (int)i, NULL, NULL);
    nodes[i]->scaler_index -= n_taxa;
    nodes[i]->next->scaler_index -= n_taxa;
    nodes[i]->next->next->scaler_index -= n_taxa;
  }

  /* set an inner node as return value */
  new_tree = nodes[n_taxa];

  /* build minimal tree with 3 tips and 1 inner node */
  pll_utree_connect_nodes(nodes[0], nodes[n_taxa],
                          PLL_TREE_DEFAULT_BRANCH_LENGTH);
  branches[n_placed_branches++] = nodes[n_taxa];
  pll_utree_connect_nodes(nodes[1], nodes[n_taxa]->next,
                          PLL_TREE_DEFAULT_BRANCH_LENGTH);
  branches[n_placed_branches++] = nodes[n_taxa]->next;
  pll_utree_connect_nodes(nodes[2], nodes[n_taxa]->next->next,
                          PLL_TREE_DEFAULT_BRANCH_LENGTH);
  branches[n_placed_branches++] = nodes[n_taxa]->next->next;

  for (i=3; i<n_taxa; ++i)
  {
    /* take tips iteratively */
    next_tip = nodes[i];
    next_inner = nodes[n_tip_nodes + i - 2];

    /* select random branch from the tree */
    rand_branch_id = (unsigned int) rand() % n_placed_branches;
    next_branch = branches[rand_branch_id];

    /* connect tip to selected branch */
    pll_utree_connect_nodes(next_branch->back, next_inner,
                            PLL_TREE_DEFAULT_BRANCH_LENGTH);
    pll_utree_connect_nodes(next_branch, next_inner->next,
                            PLL_TREE_DEFAULT_BRANCH_LENGTH);
    pll_utree_connect_nodes(next_tip, next_inner->next->next,
                            PLL_TREE_DEFAULT_BRANCH_LENGTH);

    if (pll_utree_is_tip (next_inner->back))
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
    branches[n_placed_branches++] = next_inner;
    branches[n_placed_branches++] = next_inner->next->next;
  }
  assert(n_placed_branches == max_branches);

  /* clean */
  free (nodes);
  free (branches);

  return (new_tree);
}

/* static functions */
static int utree_find_node_in_subtree(pll_utree_t * root,
                                                 pll_utree_t * node)
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

static int utree_traverse_apply(pll_utree_t * node,
                                int (*cb_pre_trav)(pll_utree_t *, void *),
                                int (*cb_post_trav)(pll_utree_t *, void *),
                                void *data)
{
  int retval = 1;

  if (pll_utree_is_tip(node))
  {
    if (cb_pre_trav)
      cb_pre_trav(node, data);
    if (cb_post_trav)
      retval = cb_post_trav(node, data);
    return retval;
  }

  if (cb_pre_trav && !cb_pre_trav(node,  data))
    return PLL_FAILURE;

  retval &= utree_traverse_apply(node->next->back,
                                 cb_pre_trav, cb_post_trav, data);

  retval &= utree_traverse_apply(node->next->next->back,
                                 cb_pre_trav, cb_post_trav, data);

  if (cb_post_trav)
    retval &= cb_post_trav(node,  data);

  return retval;
}

PLL_EXPORT int pll_utree_traverse_apply(pll_utree_t * root,
                                        int (*cb_pre_trav)(pll_utree_t *, void *),
                                        int (*cb_post_trav)(pll_utree_t *, void *),
                                        void *data)
{
  int retval = 1;

  if (pll_utree_is_tip(root)) return PLL_FAILURE;

  retval &= utree_traverse_apply(root->back, cb_pre_trav, cb_post_trav, data);
  retval &= utree_traverse_apply(root, cb_pre_trav, cb_post_trav, data);

  return retval;
}

PLL_EXPORT int pll_utree_is_tip(pll_utree_t * node)
{
  return (node->next == NULL);
}

PLL_EXPORT void pll_utree_set_length(pll_utree_t * edge,
                                            double length)
{
  edge->length = edge->back->length = length;
}

/**
 * compute the likelihood on a utree structure
 * if update_pmatrices or update_partials are set, p-matrices and CLVs are
 * updated before computing the likelihood.
 */
PLL_EXPORT double pll_utree_compute_lk(pll_partition_t * partition,
                                       pll_utree_t * tree,
                                       unsigned int * params_indices,
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

    pll_utree_traverse_apply(tree, 0, cb_update_matrices_partials, (void *) &parameters);
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



/******************************************************************************/
/* Static functions */

static int rtree_rollback_TBR(pll_tree_rollback_t * rollback_info)
{
  assert(UNIMPLEMENTED);
  return PLL_FAILURE;
}

static int rtree_rollback_SPR(pll_tree_rollback_t * rollback_info)
{
  assert(UNIMPLEMENTED);
  return PLL_FAILURE;
}

static int rtree_rollback_NNI(pll_tree_rollback_t * rollback_info)
{
  assert(UNIMPLEMENTED);
  return PLL_FAILURE;
}


static int utree_rollback_TBR(pll_tree_rollback_t * rollback_info)
{
  assert(!rollback_info->rooted);
  assert(rollback_info->rearrange_type == PLL_REARRANGE_TBR);

  pll_utree_t * p = (pll_utree_t *) rollback_info->TBR.bisect_edge;
  pll_utree_t * q = p->next->back;
  pll_utree_t * r = p->back->next->back;
  double reconn_length = rollback_info->TBR.reconn_edge.length;

  /* undo move */
    if (!pll_utree_TBR(p, &(rollback_info->TBR.reconn_edge), 0))
      return PLL_FAILURE;

  /* reset branches */
  pll_utree_set_length(p, reconn_length);
  pll_utree_set_length(q, rollback_info->TBR.bisect_left_bl);
  pll_utree_set_length(r, rollback_info->TBR.bisect_right_bl);
  pll_utree_set_length(p->next, rollback_info->TBR.reconn_parent_left_bl);
  pll_utree_set_length(p->next->next, rollback_info->TBR.reconn_parent_right_bl);
  pll_utree_set_length(p->back->next, rollback_info->TBR.reconn_child_left_bl);
  pll_utree_set_length(p->back->next->next, rollback_info->TBR.reconn_child_right_bl);

  return PLL_SUCCESS;
}

static int utree_rollback_SPR(pll_tree_rollback_t * rollback_info)
{
  assert(!rollback_info->rooted);
  assert(rollback_info->rearrange_type == PLL_REARRANGE_SPR);

  pll_utree_t * p = (pll_utree_t *) rollback_info->SPR.prune_edge;
  pll_utree_t * r = (pll_utree_t *) rollback_info->SPR.regraft_edge;
  pll_utree_t * z1 =  p->next->back;
  pll_utree_t * z2 =  r->back;

  /* undo move */
  if (!pll_utree_SPR(p, r, 0))
    return PLL_FAILURE;

  /* reset branches */
  pll_utree_set_length(z1, rollback_info->SPR.regraft_bl);
  pll_utree_set_length(p, rollback_info->SPR.prune_bl);
  pll_utree_set_length(r, rollback_info->SPR.prune_left_bl);
  pll_utree_set_length(z2, rollback_info->SPR.prune_right_bl);

  return PLL_SUCCESS;
}

static int utree_rollback_NNI(pll_tree_rollback_t * rollback_info)
{
  assert(!rollback_info->rooted);
  assert(rollback_info->rearrange_type == PLL_REARRANGE_NNI);

  pll_utree_t * p = rollback_info->NNI.edge;
  pll_utree_t * q = p->back;

  /* undo move */
  if (!pll_utree_NNI(p, rollback_info->NNI.type, 0))
      return PLL_FAILURE;

  /* reset branches */

  pll_utree_set_length(p, rollback_info->NNI.edge_bl);
  pll_utree_set_length(p->next, rollback_info->NNI.left_left_bl);
  pll_utree_set_length(p->next->next, rollback_info->NNI.left_right_bl);
  pll_utree_set_length(q->next, rollback_info->NNI.right_left_bl);
  pll_utree_set_length(q->next->next, rollback_info->NNI.right_right_bl);

  assert(UNIMPLEMENTED);
  return PLL_FAILURE;
}

/**
 * callback function for updating p-matrices and partials
 */
static int cb_update_matrices_partials(pll_utree_t * node, void *data)
{
  struct cb_params * st_data = (struct cb_params *) data;

  if (st_data->update_pmatrices)
  {
    unsigned int matrix_index = node->pmatrix_index;
    double branch_length = node->length;

    /* check integrity */
    assert(node->length == node->back->length);
    assert(node->pmatrix_index == node->back->pmatrix_index);

    pll_update_prob_matrices (st_data->partition,
                              st_data->params_indices,
                              &matrix_index,
                              &branch_length,
                              1);
  }

  if (st_data->update_partials && !pll_utree_is_tip(node))
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
