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
#ifndef PLL_TREE_H_
#define PLL_TREE_H_

#ifndef PLL_H_
#define PLL_H_
#include "pll.h"
#endif

/**
 * PLL Tree utils module
 * Prefix: pll_tree_, pll_utree_, pll_rtree_
 */
#define PLL_TREE_DEFAULT_BRANCH_LENGTH 0.1

/* NNI moves type
 *
 *   \edge       /left (next)
 *    \         /
 *     x------->
 *    /         \
 *   /           \right (next->next)
 */
#define PLL_NNI_LEFT    1
#define PLL_NNI_RIGHT   2

/* error codes (for this module, 3000-4000) ; B = 2^10+2^11*/
/* TBR errors (B + {2^2,2^1,2^0}) */
#define PLL_ERROR_TBR_LEAF_BISECTION   3073 // B + {001}
#define PLL_ERROR_TBR_OVERLAPPED_NODES 3074 // B + {010}
#define PLL_ERROR_TBR_SAME_SUBTREE     3075 // B + {011}
#define PLL_ERROR_TBR_MASK             3079 // B + {111}

/* NNI errors (B + {2^4,2^3}) */
#define PLL_ERROR_NNI_INVALID_MOVE     3080 // B + {01...}
#define PLL_ERROR_NNI_MASK             3096 // B + {11...}

/* SPR errors (B + {2^6,2^5}) */
#define PLL_ERROR_SPR_INVALID_NODE     3104 // B + {01...}
#define PLL_ERROR_SPR_MASK             3168 // B + {11...}

/* general errors (B + {2^8,2^7}) */
#define PLL_ERROR_INTERCHANGE_LEAF     3200 // B + {01...}
#define PLL_ERROR_INVALID_REARRAGE     3328 // B + {10...}

#define PLL_REARRANGE_SPR  0
#define PLL_REARRANGE_NNI  1
#define PLL_REARRANGE_TBR  2

typedef unsigned int pll_split_base_t;
typedef pll_split_base_t * pll_split_t;

typedef struct pll_tree_edge
{
  union
  {
    struct
    {
      pll_utree_t * parent;
      pll_utree_t * child;
    } utree;
    struct
    {
      pll_rtree_t * parent;
      pll_rtree_t * child;
    } rtree;
  } edge;
    double length;
} pll_tree_edge_t;

typedef struct
{
  int rearrange_type;
  int rooted;
  double  likelihood;

  union {
    struct {
      void * prune_edge;
      void * regraft_edge;
      double prune_bl;        //! length of the pruned branch
      double prune_left_bl;   //! length of the removed branch when pruning
      double prune_right_bl;  //! length of the removed branch when pruning
      double regraft_bl;      //! length of the splitted branch when regrafting
    } SPR;
    struct {
      void * edge;
      double left_left_bl;
      double left_right_bl;
      double right_left_bl;
      double right_right_bl;
      double edge_bl;
      int type;
    } NNI;
    struct {
      void * bisect_edge;
      pll_tree_edge_t reconn_edge;
      double bisect_left_bl;
      double bisect_right_bl;
      double reconn_parent_left_bl;
      double reconn_parent_right_bl;
      double reconn_child_left_bl;
      double reconn_child_right_bl;
    } TBR;
  };
} pll_tree_rollback_t;



/* Topological rearrangements */
/* functions at pll_tree.c */

PLL_EXPORT int pll_utree_TBR(pll_utree_t * b_edge,
                             pll_tree_edge_t * r_edge,
                             pll_tree_rollback_t * rollback_info);

PLL_EXPORT int pll_utree_SPR(pll_utree_t * p_edge,
                             pll_utree_t * r_edge,
                             pll_tree_rollback_t * rollback_info);

/* type = {PLL_NNI_NEXT, PLL_NNI_NEXTNEXT} */
PLL_EXPORT int pll_utree_NNI(pll_utree_t * edge,
                             int type,
                             pll_tree_rollback_t * rollback_info);

PLL_EXPORT int pll_tree_rollback(pll_tree_rollback_t * rollback_info);



/* Topological operations */
/* functions at utree_operations.c */

PLL_EXPORT int pll_utree_bisect(pll_utree_t * edge,
                                pll_utree_t ** parent_subtree,
                                pll_utree_t ** child_subtree);

PLL_EXPORT pll_tree_edge_t pll_utree_reconnect(pll_tree_edge_t * edge,
                                               pll_utree_t * pruned_edge);

PLL_EXPORT pll_utree_t * pll_utree_prune(pll_utree_t * edge);

PLL_EXPORT int pll_utree_regraft(pll_utree_t * edge,
                                 pll_utree_t * tree);

PLL_EXPORT int pll_utree_interchange(pll_utree_t * edge1,
                                     pll_utree_t * edge2);

PLL_EXPORT pll_utree_t * pll_utree_create_node(unsigned int clv_index,
                                           int scaler_index,
                                           char * label,
                                           void * data);

PLL_EXPORT int pll_utree_connect_nodes(pll_utree_t * parent,
                                       pll_utree_t * child,
                                       double length);



/* Topological search */
/* functions at utree_operations.c */

PLL_EXPORT int pll_utree_nodes_at_edge_dist(pll_utree_t * edge,
                                            pll_utree_t ** outbuffer,
                                            unsigned int * n_nodes,
                                            unsigned int distance,
                                            int fixed);


PLL_EXPORT int pll_utree_nodes_at_node_dist(pll_utree_t * node,
                                            pll_utree_t ** outbuffer,
                                            unsigned int * n_nodes,
                                            unsigned int distance,
                                            int fixed);



/* Tree construction */
/* functions at pll_tree.c */

PLL_EXPORT pll_utree_t * pll_utree_create_random(unsigned int n_taxa,
                                                 const char * const* names);



/* Discrete operations */
/* functions at utree_distances.c */

PLL_EXPORT unsigned int pll_utree_rf_distance(pll_utree_t * t1,
                                              pll_utree_t * t2,
                                              unsigned int n_tips);

PLL_EXPORT unsigned int pll_utree_rf_split_distance(pll_split_t * s1,
                                                    pll_split_t * s2,
                                                    unsigned int n_tips);

PLL_EXPORT pll_split_t * pll_utree_split_create(pll_utree_t * tree,
                                                unsigned int n_tips,
                                                unsigned int * n_splits);

PLL_EXPORT void pll_utree_show_split(pll_split_t split, unsigned int n_tips);

PLL_EXPORT void pll_utree_split_destroy(pll_split_t * split_list);



/* Additional utilities */
/* functions at pll_tree.c */

PLL_EXPORT int pll_utree_traverse_apply(pll_utree_t * root,
                                        int (*cb_pre_trav)(pll_utree_t *, void *),
                                        int (*cb_post_trav)(pll_utree_t *, void *),
                                        void *data);

PLL_EXPORT inline int pll_utree_is_tip(pll_utree_t * node);

PLL_EXPORT inline void pll_utree_set_length(pll_utree_t * edge,
                                            double length);

PLL_EXPORT double pll_utree_compute_lk(pll_partition_t * partition,
                                       pll_utree_t * tree,
                                       unsigned int * params_indices,
                                       int update_pmatrices,
                                       int update_partials);

#endif /* PLL_TREE_H_ */
