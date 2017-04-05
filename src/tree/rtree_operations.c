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

/**
 * @file rtree_operations.c
 *
 * @brief Operations on rooted tree structures
 *
 * @author Diego Darriba
 */

/* finds the self and sister pointers */
PLL_EXPORT int pllmod_rtree_get_sibling_pointers(pll_rnode_t * node,
                                                 pll_rnode_t ***self,
                                                 pll_rnode_t ***sister)
{
  if (!node->parent)
  {
    /* if there is no parent, there are no pointers */
    if (self) *self = NULL;
    if (sister) *sister = NULL;
  }
  else if (node->parent->left == node)
  {
    if (self) *self = &(node->parent->left);
    if (sister) *sister = &(node->parent->right);
  }
  else if (node->parent->right == node)
  {
    if (self) *self = &(node->parent->right);
    if (sister) *sister = &(node->parent->left);
  }
  else
  {
    /* `node` is not the left nor the right child of its parent */
    if (self) *self = NULL;
    if (sister) *sister = NULL;
    pllmod_set_error(PLLMOD_TREE_ERROR_INVALID_TREE,
                     "Tree is not consistent");
    return PLL_FAILURE;
  }

  return PLL_SUCCESS;
}

/**
 * @brief Prunes a subtree in a rooted tree
 * @param node the node to prune
 * @return the new connected node, if the operation was applied correctly, NULL otherwise
 */
PLL_EXPORT pll_rnode_t * pllmod_rtree_prune(pll_rnode_t * node)
{
  pll_rnode_t **self_ptr, **sister_ptr, **parent_ptr;
  pll_rnode_t *connected_node = NULL;
  assert(node);

  if (!node->parent)
  {
    pllmod_set_error(PLLMOD_TREE_ERROR_SPR_INVALID_NODE,
                     "Attempting to prune the root node");
    return NULL;
  }

  if (!pllmod_rtree_get_sibling_pointers(node,
                               &self_ptr,
                               &sister_ptr))
  {
    /* return and spread error */
    return NULL;
  }
  else
  {
    assert (self_ptr && sister_ptr);
  }

  /* connect adjacent subtrees together */
  if (node->parent->parent)
  {
    /* connect parent->parent and sister */
    connected_node = node->parent->parent;
    if (!pllmod_rtree_get_sibling_pointers(node->parent,
                                 &parent_ptr,
                                 NULL))
    {
      /* return and spread error */
      return NULL;
    }
    else
    {
      assert (parent_ptr);
    }
    *parent_ptr = *sister_ptr;
    (*sister_ptr)->parent = node->parent->parent;

    /* disconnect pruned tree */
    *sister_ptr = NULL;
    node->parent->parent = NULL;
  }
  else
  {
    /* re-root */
    connected_node = *sister_ptr;

    (*sister_ptr)->parent = NULL;

    /* disconnect pruned tree */
    *sister_ptr = NULL;
  }

  return connected_node;
}

/**
 * @brief Regrafts a dettached subtree into a branch
 * 
 * @param node the node to regraft
 * @param tree the target branch
 *
 * @return PLL_SUCCESS if the operation was applied correctly,
 *         PLL_FAILURE otherwise (check pll_errmsg for details)
 */
PLL_EXPORT int pllmod_rtree_regraft(pll_rnode_t * node,
                                    pll_rnode_t * tree)
{
  pll_rnode_t *parent_node;
  pll_rnode_t **edge_from_parent = 0,
              **edge_to_child    = 0;

  /* node must contain a dettached parent */
  if (!node->parent || node->parent->parent)
  {
    pllmod_set_error(PLLMOD_TREE_ERROR_SPR_INVALID_NODE,
                     "Attempting to regraft a node without dettached parent");
    return PLL_FAILURE;
  }

  /* `node` parent should contain a pointer to `node` */
  if (!pllmod_rtree_get_sibling_pointers(node,
                               NULL,
                               &edge_to_child))
  {
    /* return and spread error */
    assert(pll_errno);
    return PLL_FAILURE;
  }

  parent_node = node->parent;

  if (tree->parent)
  {
    if (tree->parent && !pllmod_rtree_get_sibling_pointers(tree,
                                  &edge_from_parent,
                                  NULL))
    {
      /* return and spread error */
      assert(pll_errno);
      return PLL_FAILURE;
    }
    assert(edge_from_parent);

    /* set new parents */
    parent_node->parent = tree->parent;

    /* set new children */
    *edge_from_parent = parent_node;
  }

  tree->parent = parent_node;
  *edge_to_child = tree;

  return PLL_SUCCESS;
}

/**
 * Performs one SPR move
 * The CLV, scaler and pmatrix indices are updated.
 *
 * @param[in] p_node Edge to be pruned
 * @param[in] r_tree Edge to be regrafted
 * @param[in] root The tree root (it might change)
 * @param[in,out] rollback_info Rollback information
 * @return PLL_SUCCESS if the move was applied correctly,
 *         PLL_FAILURE otherwise (check pll_errmsg for details)
 */
PLL_EXPORT int pllmod_rtree_spr(pll_rnode_t * p_node,
                                pll_rnode_t * r_tree,
                                pll_rnode_t ** root,
                                pll_tree_rollback_t * rollback_info)
{
  pll_rnode_t **self_ptr, **sister_ptr;

  if (!pllmod_rtree_get_sibling_pointers(p_node,
                               &self_ptr,
                               &sister_ptr))
  {
    /* return and spread error */
    assert(pll_errno);
    return PLL_FAILURE;
  }
  else
  {
    assert (self_ptr && sister_ptr);
  }

  /* save rollback information */
  if (rollback_info)
  {
    rollback_info->rearrange_type     = PLLMOD_TREE_REARRANGE_SPR;
    rollback_info->rooted             = 1;
    rollback_info->SPR.prune_edge     = (void *) p_node;
    rollback_info->SPR.regraft_edge   = (void *) *sister_ptr;
    //TODO: Set branch lengths
    // rollback_info->SPR.prune_bl       = p_edge->parent->length;
    // rollback_info->SPR.prune_left_bl  = (*sister)->length;
    // rollback_info->SPR.prune_right_bl = p_edge->->length;
    // rollback_info->SPR.regraft_bl     = r_tree->length;
  }

  if (pllmod_rtree_prune(p_node) == NULL)
  {
    /* return and spread error */
    assert(pll_errno);
    return PLL_FAILURE;
  }

  if (!pllmod_rtree_regraft(p_node,
                            r_tree))
  {
    /* return and spread error */
    assert(pll_errno);
    return PLL_FAILURE;
  }

  /* reset root in case it has changed */
  while (root && (*root)->parent) *root = (*root)->parent;

  return PLL_SUCCESS;
}

/* re-roots the tree at the branch connecting `new_root` with its parent */
PLL_EXPORT int pllmod_rtree_reroot(pll_rnode_t * root,
                                   pll_rnode_t * new_root)
{
  return PLL_FAILURE;
}

static void rtree_nodes_at_node_dist_down(pll_rnode_t * root,
                                          pll_rnode_t ** outbuffer,
                                          unsigned int * node_count,
                                          int min_distance,
                                          int max_distance)
{
  if (max_distance < 0) return;

  if (min_distance < 0)
  {
    outbuffer[*node_count] = root;
    *node_count = *node_count + 1;
  }

  if (!(root->left && root->right)) return;

  rtree_nodes_at_node_dist_down(root->left,
                                outbuffer,
                                node_count,
                                min_distance-1,
                                max_distance-1);
  rtree_nodes_at_node_dist_down(root->right,
                                outbuffer,
                                node_count,
                                min_distance-1,
                                max_distance-1);
}

PLL_EXPORT int pllmod_rtree_nodes_at_node_dist(pll_rnode_t * root,
                                               pll_rnode_t ** outbuffer,
                                               unsigned int * node_count,
                                               int min_distance,
                                               int max_distance)
{
  pll_rnode_t * current_root = root;
  pll_rnode_t ** sister_ptr;

  if (max_distance < min_distance)
    {
      pllmod_set_error(PLLMOD_ERROR_INVALID_RANGE,
                 "Invalid distance range: %d..%d (max_distance < min_distance)",
                 min_distance, max_distance);
      return PLL_FAILURE;
    }

  *node_count = 0;

  while (current_root->parent)
  {
    if (!pllmod_rtree_get_sibling_pointers(current_root,
                                 NULL,
                                 &sister_ptr))
    {
      /* return and spread error */
      assert(pll_errno);
      return PLL_FAILURE;
    }

    --min_distance;
    --max_distance;

    current_root = current_root->parent;
    if (min_distance < 0)
    {
      outbuffer[*node_count] = current_root;
      *node_count = *node_count + 1;
    }
    rtree_nodes_at_node_dist_down(*sister_ptr,
                                  outbuffer,
                                  node_count,
                                  min_distance-1,
                                  max_distance-1);
  }


  return PLL_SUCCESS;
}
