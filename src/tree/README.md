# libpll/tree

[![License](https://img.shields.io/badge/license-AGPL-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

## Introduction

Module for tree management. Exportable functions in this module start either with the prefix `pllmod_tree_` for general tree management functions, `pllmod_utree_` for unrooted tree operations, or `pllmod_rtree_` for rooted tree operations.

## Compilation instructions

Read the compilation instructions at the main README file

## Code

The code is written in C.

|    File               | Description                                   |
|-----------------------|-----------------------------------------------|
|**pll_tree.c**         | Functions for performing complete operations. |
|**utree_operations.c** | Operations on unrooted trees.                 |
|**rtree_operations.c** | Operations on rooted trees.                   |
|**tree_hashtable.c**   | Operations on unrooted trees.                 |
|**consensus.c**        | Functions for consensus trees.                |
|**treeinfo.c**         | Functions related to global tree information. |

## Type definitions

* unsigned int `pll_split_base_t`
* `pll_split_base_t * pll_split_t`
* struct `pll_split_system_t`
* struct `pll_tree_rollback_t`
* struct `pllmod_treeinfo_t`

## Flags

* `PLLMOD_TREE_REARRANGE_SPR`
* `PLLMOD_TREE_REARRANGE_NNI`
* `PLLMOD_TREE_REARRANGE_TBR`

* `PLLMOD_TREE_BRLEN_LINKED`
* `PLLMOD_TREE_BRLEN_SCALED`
* `PLLMOD_TREE_BRLEN_UNLINKED`

* `PLLMOD_TREEINFO_PARTITION_ALL`

## Functions

* `int pllmod_utree_tbr`
* `int pllmod_utree_spr`
* `int pllmod_utree_nni`
* `int pllmod_tree_rollback`
* `int pllmod_rtree_spr`
* `int pllmod_rtree_get_sibling_pointers`
* `pll_rtree_t * pllmod_rtree_prune`
* `int pllmod_rtree_regraft`
* `int pllmod_utree_bisect`
* `pll_tree_edge_t pllmod_utree_reconnect`
* `pll_utree_t * pllmod_utree_prune`
* `int pllmod_utree_regraft`
* `int pllmod_utree_interchange`
* `pll_utree_t * pllmod_utree_create_node`
* `int pllmod_utree_connect_nodes`
* `int pllmod_rtree_nodes_at_node_dist`
* `int pllmod_utree_nodes_at_node_dist`
* `int pllmod_utree_nodes_at_edge_dist`
* `pll_utree_t * pllmod_utree_create_random`
* `unsigned int pllmod_utree_rf_distance`
* `int pllmod_utree_consistency_check`
* `int pllmod_utree_consistency_set`
* `unsigned int pllmod_utree_split_rf_distance`
* `pll_split_t * pllmod_utree_split_create`
* `void pllmod_utree_split_normalize_and_sort`
* `void pllmod_utree_split_show`
* `void pllmod_utree_split_destroy`
* `int pllmod_utree_compatible_splits`
* `pll_utree_t * pllmod_utree_from_splits`
* `pll_utree_t * pllmod_utree_consensus`
* `int pllmod_utree_set_clv_minimal`
* `int pllmod_utree_traverse_apply`
* `int pllmod_utree_is_tip`
* `void pllmod_utree_set_length`
* `void pllmod_utree_scale_branches`
* `double pllmod_utree_compute_lk`
* `int pllmod_rtree_traverse_apply`
* `pllmod_treeinfo_t * pllmod_treeinfo_create`
* `int pllmod_treeinfo_init_partition`
* `int pllmod_treeinfo_set_active_partition`
* `void pllmod_treeinfo_set_root`
* `void pllmod_treeinfo_set_branch_length`
* `int pllmod_treeinfo_destroy_partition`
* `void pllmod_treeinfo_destroy`
* `int pllmod_treeinfo_update_prob_matrices`
* `void pllmod_treeinfo_invalidate_all`
* `int pllmod_treeinfo_validate_clvs`
* `void pllmod_treeinfo_invalidate_pmatrix`
* `void pllmod_treeinfo_invalidate_clv`
* `double pllmod_treeinfo_compute_loglh`

## Error codes

* 3073: `PLLMOD_TREE_ERROR_TBR_LEAF_BISECTION`
* 3074: `PLLMOD_TREE_ERROR_TBR_OVERLAPPED_NODES`
* 3075: `PLLMOD_TREE_ERROR_TBR_SAME_SUBTREE`
* 3079: `PLLMOD_TREE_ERROR_TBR_MASK`

* 3080: `PLLMOD_TREE_ERROR_NNI_INVALID_MOVE`
* 3096: `PLLMOD_TREE_ERROR_NNI_MASK`

* 3104: `PLLMOD_TREE_ERROR_SPR_INVALID_NODE`
* 3168: `PLLMOD_TREE_ERROR_SPR_MASK`

* 3200: `PLLMOD_TREE_ERROR_INTERCHANGE_LEAF`
* 3328: `PLLMOD_TREE_ERROR_INVALID_REARRAGE`
* 3456: `PLLMOD_TREE_ERROR_INVALID_TREE_SIZE`
* 3584: `PLLMOD_TREE_ERROR_INVALID_TREE`
* 3712: `PLLMOD_TREE_ERROR_INVALID_SPLIT`
* 3840: `PLLMOD_TREE_ERROR_EMPTY_SPLIT`
* 3968: `PLLMOD_TREE_ERROR_INVALID_THRESHOLD`
