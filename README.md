# pll-modules
High Level modules for the Low Level Phylogenetic Likelihood Library

## Introduction

## Compilation instructions

PLL-Modules depends on the PLL library submodule (libs/libpll)

## Usage examples 

Please refer to the [wiki page](https://github.com/ddarriba/pll-modules/wiki) and/or the [examples directory](https://github.com/ddarriba/pll-modules/tree/master/examples).

## Documentation

Please refer to the [wiki page](https://github.com/ddarriba/pll-modules/wiki).

## Available functionality

Below is a list of available functions in the current version.

### Tree utilities module (src/tree)

This module contains functions for manipulating a tree structure (utree/rtee).
Prefixes: `pll_tree_` (general), `pll_utree_` (unrooted tree) or `pll_rtree_` (rooted tree). 

* `int pll_utree_TBR(pll_utree_t * b_edge, pll_tree_edge_t * r_edge);`
* `pll_utree_t * pll_utree_create_random(unsigned int n_taxa, const char ** names);`
* `pll_utree_t * pll_utree_create_node(unsigned int clv_index, int scaler_index, char * label, void * data);`
* `int pll_utree_connect_nodes(pll_utree_t * parent, pll_utree_t * child, double length);`
* `int pll_utree_bisect(pll_utree_t * edge, pll_tree_edge_t * parent_subtree, pll_tree_edge_t * child_subtree);`
* `pll_tree_edge_t pll_utree_reconnect(pll_tree_edge_t * edge, unsigned int parent_pmatrix_index, unsigned int parent_clv_index, int parent_scaler_index, unsigned int child_pmatrix_index, unsigned int child_clv_index, int child_scaler_index, unsigned int edge_pmatrix_index);`
* `int pll_utree_nodes_at_edge_dist(pll_utree_t * edge, pll_utree_t ** outbuffer, unsigned int * n_nodes, unsigned int distance, int fixed);`
* `int pll_utree_nodes_at_node_dist(pll_utree_t * node, pll_utree_t ** outbuffer, unsigned int * n_nodes, unsigned int distance, int fixed);`
* `int pll_utree_interchange(pll_utree_t * edge1, pll_utree_t * edge2);`

### MSA utilities module (src/msa)

This module contains functions for analyzing MSAs.
Prefixes: `pll_msa_`

* `double * pll_msa_empirical_frequencies(pll_partition_t * partition);`
* `double * pll_msa_empirical_subst_rates(pll_partition_t * partition);`
* `double pll_msa_empirical_invariant_sites(pll_partition_t *partition);`

### Optimization module (src/optimize)

This module contains functions for model parameter optimization (libpll_optimize).
Prefixes: `pll_optimize_` (high level) or `pll_minimize_` (low level).

#### High level optimization functions

* `double pll_optimize_parameters_onedim(pll_optimize_options_t * p, double min, double max);`
* `double pll_optimize_parameters_multidim(pll_optimize_options_t * p, double *umin, double *umax);`
* `double pll_optimize_branch_lengths_iterative (pll_partition_t * partition, pll_utree_t * tree, unsigned int params_index, unsigned int freqs_index, double tolerance, int smoothings, int keep_update);`
* `double pll_optimize_branch_lengths_local (pll_partition_t * partition, pll_utree_t * tree, unsigned int params_index, unsigned int freqs_index, double tolerance, int smoothings, int radius, int keep_update);`

#### Low level optimization functions

* `double pll_minimize_newton(double x1, double xguess, double x2, unsigned int max_iters, double *score, void *params, double (deriv_func)(void *, double, double *, double *));`
* `double pll_minimize_lbfgsb(double *x, double *xmin, double *xmax, int *bound, unsigned int n, double factr, double pgtol, void *params, double (*target_funk)(void *, double *));`
* `double pll_minimize_brent(double xmin, double xguess, double xmax, double xtol, double *fx, double *f2x, void *params, double (*target_funk)(void *, double));`
