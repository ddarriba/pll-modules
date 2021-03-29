# libpll/optimize

[![License](https://img.shields.io/badge/license-AGPL-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

## Introduction

Module for model parameters optimization. Exportable functions in this module start either with the prefix `pllmod_optimize_` for high level optimization algorithms, or `pllmod_minimize_` for low level optimization functions.

## Compilation instructions

Read the compilation instructions at the main README file

## Code

The code is written in C.

|    File              | Description                                |
|----------------------|--------------------------------------------|
|**pll_optimize.c**    | High level optimization algorithms.        |
|**opt_algorithms.c**  | Low level optimization algorithms.         |
|**lbfgsb/**           | Core implementation of L-BFGS-B algorithm. |

## Type definitions

* struct `pll_likelihood_info_t`
* struct `pll_optimize_options_t`
* struct `pll_newton_tree_params_t`
* struct `pll_newton_tree_params_multi_t`

## Flags

* `PLLMOD_OPT_PARAM_SUBST_RATES`
* `PLLMOD_OPT_PARAM_ALPHA`
* `PLLMOD_OPT_PARAM_PINV`
* `PLLMOD_OPT_PARAM_FREQUENCIES`
* `PLLMOD_OPT_PARAM_BRANCHES_SINGLE`
* `PLLMOD_OPT_PARAM_BRANCHES_ALL`
* `PLLMOD_OPT_PARAM_BRANCHES_ITERATIVE`
* `PLLMOD_OPT_PARAM_TOPOLOGY`
* `PLLMOD_OPT_PARAM_FREE_RATES`
* `PLLMOD_OPT_PARAM_RATE_WEIGHTS`
* `PLLMOD_OPT_PARAM_BRANCH_LEN_SCALER`
* `PLLMOD_OPT_PARAM_ALL`

## Functions

* `double pllmod_opt_minimize_newton`
* `double pllmod_opt_minimize_lbfgsb`
* `double pllmod_opt_minimize_brent`
* `void pllmod_opt_minimize_em`
* `void pllmod_opt_derivative_func`
* `double pllmod_opt_optimize_onedim`
* `double pllmod_opt_optimize_multidim`
* `double pllmod_opt_optimize_branch_lengths_iterative`
* `double pllmod_opt_optimize_branch_lengths_local`
* `double pllmod_opt_optimize_branch_lengths_local_multi`

## Error codes

* 2000: `PLLMOD_OPT_ERROR_PARAMETER`
* 2010: `PLLMOD_OPT_ERROR_TAXA_MISMATCH`
* 2020: `PLLMOD_OPT_ERROR_SEQLEN_MISMATCH`
* 2030: `PLLMOD_OPT_ERROR_ALIGN_UNREADABLE`
* 2100: `PLLMOD_OPT_ERROR_LBFGSB_UNKNOWN`
* 2210: `PLLMOD_OPT_ERROR_NEWTON_DERIV`
* 2220: `PLLMOD_OPT_ERROR_NEWTON_LIMIT`
* 2230: `PLLMOD_OPT_ERROR_NEWTON_UNKNOWN`
