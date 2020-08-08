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
  * @file rtree_operations.c
  *
  * @brief Operations on unrooted tree structures
  *
  * @author Diego Darriba
  */

#include "pll_tree.h"

#include "../pllmod_common.h"

static void utree_nodes_at_dist(pll_unode_t * node,
                                pll_unode_t ** outbuffer,
                                unsigned int * index,
                                unsigned int min_distance,
                                unsigned int max_distance,
                                unsigned int depth);



/******************************************************************************/
/* Topological operations */

/**
 * @brief Bisects the tree by removing one edge
 *
 * Removes the edge \p edge and frees the nodes defining that edge.
 * Reconnects the subtrees at the sides of the edge (figure below).
 * The branch lengths of the new edges are the sum of the removed ones.
 * The join branch contains the pmatrix index of the parent edges
 * The removed pmatrix indices are returned in the field
 *     'additional_pmatrix_index' of both output subtrees
 *
 * Returns the new parent and child edges, where parent is the closest to \p edge.
 *
 *   A            C              A        C
 *    \___edge___/       ---->   |        |
 *    /          \               |        |
 *   B            D              B        D
 *   A,B,C,D are subtrees
 *
 * @param[in] edge            edge to remove
 * @param[out] parent_subtree edge corresponding to the 'edge' subtree
 * @param[out] child_subtree  edge corresponding to the 'edge->back' subtree
 * @return PLL_SUCCESS if OK
 */
PLL_EXPORT int pllmod_utree_bisect(pll_unode_t * edge,
                                   pll_unode_t ** parent_subtree,
                                   pll_unode_t ** child_subtree)
{
  assert(parent_subtree);
  assert(child_subtree);

  pll_unode_t * aux_tree;

  if (!edge->next)
    return PLL_FAILURE;

  pll_unode_t * c_edge = edge->back;

  /* connect parent subtree */
  (*parent_subtree) = edge->next->back;
  aux_tree = edge->next->next->back;

  pllmod_utree_connect_nodes(*parent_subtree,
                             aux_tree,
                             (*parent_subtree)->length + aux_tree->length);

  edge->next->pmatrix_index = edge->next->next->pmatrix_index;

  /* connect child subtree */
  (*child_subtree) = c_edge->next->back;
  aux_tree = c_edge->next->next->back;

  pllmod_utree_connect_nodes(*child_subtree,
                             aux_tree,
                             (*child_subtree)->length + aux_tree->length);

  c_edge->next->pmatrix_index = c_edge->next->next->pmatrix_index;

  return PLL_SUCCESS;
}

/**
 * Reconnects two subtrees by adding 2 new nodes and 1 edge.
 *
 * Adds 1 new edge connecting edges \p edge.parent and \p edge.child with
 * length \p edge.length.
 *
 *   A       C         A              C
 *   |       |  ---->   \            /
 *                       e1--edge--e2
 *   |       |          /            \
 *   B       D         B              D
 *   A,B,C,D are subtrees
 *
 * @param edge                 new edge (edge structure)
 * @param pruned_edge          edge to prune, defined by a tree node
 *
 * @return the new created edge
 */
PLL_EXPORT pll_tree_edge_t pllmod_utree_reconnect(pll_tree_edge_t * edge,
                                                  pll_unode_t * pruned_edge)
{
  /* create and connect 2 new nodes */
  pll_unode_t *parent_node, *child_node;
  assert(pruned_edge->back);

  parent_node = pruned_edge;
  child_node  = pruned_edge->back;
  assert(parent_node->back == child_node && child_node->back == parent_node);

  assert(!pllmod_utree_is_tip(parent_node));
  assert(!pllmod_utree_is_tip(child_node));

  pll_tree_edge_t new_edge;
  new_edge.edge.utree.child = child_node;
  new_edge.length = edge->length;

  /* set length */
  pllmod_utree_set_length(parent_node, edge->length);

  /* reconnect parent close to edge.parent */
  pllmod_utree_connect_nodes(parent_node->next->next,
                             edge->edge.utree.parent->back,
                             edge->edge.utree.parent->back->length);

  pllmod_utree_connect_nodes(edge->edge.utree.parent,
                             parent_node->next,
                             0);

  /* reconnect child close to edge.child */
  pllmod_utree_connect_nodes(child_node->next->next,
                             edge->edge.utree.child->back,
                             edge->edge.utree.child->back->length);

  pllmod_utree_connect_nodes(edge->edge.utree.child,
                             child_node->next,
                             0);

  return new_edge;
}

/**
 * @brief Prunes a subtree in an unrooted tree
 *
 * Disconnecs an edge (e1) and connects the adjacent nodes. New branch (A-B)
 * length is set to the sum of lengths of previous branch (e1-A + e1-B)
 *
 *   A              C              A                   C
 *    \            /               |                  /
 *     e1--edge--e2        --->    |  +   e1--edge--e2
 *    /            \               |                  \
 *   B              D              B                   D
 *   A,B,C,D are subtrees
 *
 *  Note that `edge` is disconnected after the operation
 *
 * @param edge the edge to prune
 * @return the new connected edge, if the operation was applied correctly
 */
PLL_EXPORT pll_unode_t * pllmod_utree_prune(pll_unode_t * edge)
{
  pll_unode_t *edge1, *edge2;

  assert(edge);
  if (!edge->next)
  {
    /* invalid node */
    pllmod_set_error(PLLMOD_TREE_ERROR_SPR_INVALID_NODE,
                     "Attempting to prune a tip node");
    return NULL;
  }

  /* connect adjacent subtrees together */
  edge1 = edge->next->back;
  edge2 = edge->next->next->back;
  pllmod_utree_connect_nodes(edge1, edge2, edge1->length + edge2->length);

  /* disconnect pruned edge */
  edge->next->back = edge->next->next->back = NULL;

  return edge1;
}

/**
 * @brief Regrafts an edge into a tree
 *
 * Connects a disconnected edge (provided by `e2` in the graph below)
 * into a tree
 *
 *  A                    C         A              C
 *   \                   |          \            /
 *    e1--edge--e2   +   |   --->    e1--edge--e2
 *   /                   |          /            \
 *  B                    D         B              D
 *   A,B,C,D are subtrees
 *
 *  The length of the new branches (e2-C and e2-D) are set to half the length
 *  of the removed branch (C-D)
 *
 * @param edge the edge to regraft
 * @param tree the tree to connect `edge` to
 * @return PLL_SUCCESS if the operation was applied correctly,
 *         PLL_FAILURE otherwise (check pll_errmsg for details)
 */
PLL_EXPORT int pllmod_utree_regraft(pll_unode_t * edge,
                                    pll_unode_t * tree)
{
  pll_unode_t *edge1, *edge2;
  double new_length;

  assert(edge && tree);
  if (!edge->next)
  {
    /* invalid node */
    pllmod_set_error(PLLMOD_TREE_ERROR_SPR_INVALID_NODE,
                     "Attempting to regraft a tip node");
    return PLL_FAILURE;
  }
  if (edge->next->back || edge->next->next->back)
  {
    /* invalid node */
    pllmod_set_error(PLLMOD_TREE_ERROR_SPR_INVALID_NODE,
                     "Attempting to regraft a connected node");
    return PLL_FAILURE;
  }

  /* connect tree with edge, splitting the branch designed by tree */
  edge1      = tree;
  edge2      = tree->back;
  new_length = tree->length/2;
  pllmod_utree_connect_nodes(edge1, edge->next,       new_length);
  pllmod_utree_connect_nodes(edge->next->next, edge2, new_length);

  return PLL_SUCCESS;
}

/**
 * @brief Interchanges 2 edges, represented by 2 internal nodes
 *
 * CLV and scaler indices, and labels are interchanged between nodes to match
 * the other 2 nodes in the triplet.
 *
 * @return PLL_SUCCESS if the operation was applied correctly,
 *         PLL_FAILURE otherwise (check pll_errmsg for details)
 */
PLL_EXPORT int pllmod_utree_interchange(pll_unode_t * node1,
                                        pll_unode_t * node2)
{
  pll_unode_t *next1 = node2->back;
  pll_unode_t *next2 = node1->back;

  pllmod_utree_connect_nodes(next1, node1, next1->length);
  pllmod_utree_connect_nodes(next2, node2, next2->length);

  return PLL_SUCCESS;
}

/**
 * @brief Creates a new circular node
 *
 *           n2
 *          / |
 *        n1  |
 *          \ |
 *           n3
 *
 * All parameters are shared among the nodes in the triplet
 *
 * @param clv_index    the clv_index
 * @param scaler_index the scaler index
 * @param label        the node label
 * @param data         the data pointer
 *
 * @return the new node
 */
PLL_EXPORT pll_unode_t * pllmod_utree_create_node(unsigned int clv_index,
                                                  int scaler_index,
                                                  char * label,
                                                  void * data)
{
  pll_unode_t * new_node = (pll_unode_t *)calloc(1, sizeof(pll_unode_t));
  new_node->next         = (pll_unode_t *)calloc(1, sizeof(pll_unode_t));
  new_node->next->next   = (pll_unode_t *)calloc(1, sizeof(pll_unode_t));
  if (!(new_node && new_node->next && new_node->next->next))
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for new node\n");
    return NULL;
  }

  new_node->next->next->next = new_node;
  new_node->label = label;
  new_node->next->label =
    new_node->next->next->label =
    new_node->label;
  new_node->next->data =
    new_node->next->next->data =
    new_node->data = data;
  new_node->next->length =
    new_node->next->next->length =
    new_node->length = 0;
  new_node->next->clv_index =
    new_node->next->next->clv_index =
    new_node->clv_index = clv_index;
  new_node->next->scaler_index =
    new_node->next->next->scaler_index =
    new_node->scaler_index = scaler_index;
  new_node->back =
    new_node->next->back =
    new_node->next->next->back = NULL;
  return new_node;
}

/**
 * @brief Connects 2 nodes and sets the pmatrix index and branch length
 *
 * connects `back` pointers of `child` and `parent`
 * pmatrix index for `child` is set to the one in `parent`
 *
 * @param[in,out] parent the parent node
 * @param[in,out] child  the child node
 * @param[in] length     the branch length
 *
 */
PLL_EXPORT int pllmod_utree_connect_nodes(pll_unode_t * parent,
                                          pll_unode_t * child,
                                          double length)
{
  if(!(parent && child))
    return PLL_FAILURE;

  parent->back = child;
  child->back = parent;
  pllmod_utree_set_length(parent, length);

  /* PMatrix index is set to parent node */
  child->pmatrix_index = parent->pmatrix_index;

  return PLL_SUCCESS;
}

/******************************************************************************/
/* Topological search */

/**
 * Returns the list of nodes at a distance between \p min_distance and
 * \p max_distance from a specified node
 *
 * @param[in] node the root node
 * @param[out] outbuffer the list of nodes. Outbuffer should be allocated
 * @param[out] node_count the number of nodes returned in \p outbuffer
 * @param[in] min_distance the minimum distance to check
 * @param[in] max_distance the maximum distance to check
 */
PLL_EXPORT int pllmod_utree_nodes_at_node_dist(pll_unode_t * node,
                                               pll_unode_t ** outbuffer,
                                               unsigned int * node_count,
                                               unsigned int min_distance,
                                               unsigned int max_distance)
{
  if (!node->next)
  {
    pllmod_set_error(PLLMOD_ERROR_INVALID_NODE_TYPE,
                     "Internal node expected, but tip node was provided");
    return PLL_FAILURE;
  }

  if (max_distance < min_distance)
    {
      pllmod_set_error(PLLMOD_ERROR_INVALID_RANGE,
                 "Invalid distance range: %d..%d (max_distance < min_distance)",
                 min_distance, max_distance);
      return PLL_FAILURE;
    }

  *node_count = 0;

  /* we will traverse an unrooted tree in the following way

               1
             /
          --*
             \
               2
    */


  utree_nodes_at_dist(node, outbuffer, node_count, min_distance, max_distance, 0);

  return PLL_SUCCESS;
}

/**
 * Returns the list of nodes at a distance between \p min_distance and
 * \p max_distance from a specified edge
 *
 * @param[in] edge the root edge
 * @param[out] outbuffer the list of nodes. Outbuffer should be allocated
 * @param[out] node_count the number of nodes returned in \p outbuffer
 * @param[in] min_distance the minimum distance to check
 * @param[in] max_distance the maximum distance to check
 */

PLL_EXPORT int pllmod_utree_nodes_at_edge_dist(pll_unode_t * edge,
                                               pll_unode_t ** outbuffer,
                                               unsigned int * node_count,
                                               unsigned int min_distance,
                                               unsigned int max_distance)
{
  unsigned int depth = 0;

  if (!edge->next)
  {
    pllmod_set_error(PLLMOD_ERROR_INVALID_NODE_TYPE,
                     "Internal node expected, but tip node was provided");
    return PLL_FAILURE;
  }

  if (max_distance < min_distance)
    {
      pllmod_set_error(PLLMOD_ERROR_INVALID_RANGE,
                 "Invalid distance range: %d..%d (max_distance < min_distance)",
                 min_distance, max_distance);
      return PLL_FAILURE;
    }

  *node_count = 0;

  /* we will traverse an unrooted tree in the following way

       3          1
        \        /
         * ---- *
        /        \
       4          2
   */

  utree_nodes_at_dist(edge->back, outbuffer, node_count,
                      min_distance, max_distance, depth+1);
  utree_nodes_at_dist(edge, outbuffer, node_count,
                      min_distance, max_distance, depth);

  return PLL_SUCCESS;
}


/******************************************************************************/
/* static functions */

static void utree_nodes_at_dist(pll_unode_t * node,
                                pll_unode_t ** outbuffer,
                                unsigned int * index,
                                unsigned int min_distance,
                                unsigned int max_distance,
                                unsigned int depth)
{
  if (depth >= min_distance && depth <= max_distance)
  {
    outbuffer[*index] = node;
    *index = *index + 1;
  }

  if (depth >= max_distance || !(node->next)) return;

  utree_nodes_at_dist(node->next->back, outbuffer, index,
                      min_distance, max_distance, depth+1);
  utree_nodes_at_dist(node->next->next->back, outbuffer, index,
                      min_distance, max_distance, depth+1);
}
