# libpll/algorithm

[![License](https://img.shields.io/badge/license-AGPL-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

## Introduction

Module for miscellaneous algorithms involving several modules together
Exportable functions in this module start either with the prefix `pllmod_algo_`.

## Compilation instructions

Read the compilation instructions at the main README file

## Code

The code is written in C.

|     File              | Description                                |
|-----------------------|--------------------------------------------|
|**pllmod_algorithm.c** | High level algorithms.                     |
|**algo_callback.c**    | Internal callback functions.               |
|**algo_search.c**      | Internal functions for topological search. |

## Type definitions

* struct `cutoff_info_t`

## Functions

### Functions for optimization

* `double pllmod_algo_opt_frequencies`
* `double pllmod_algo_opt_subst_rates`
* `double pllmod_algo_opt_alpha`
* `double pllmod_algo_opt_pinv`
* `double pllmod_algo_opt_alpha_pinv`
* `double pllmod_algo_opt_rates_weights`
* `double pllmod_algo_opt_brlen_scaler`

### Functions for topological search

* `double pllmod_algo_spr_round`
