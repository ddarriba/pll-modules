# libpll/msa

[![License](https://img.shields.io/badge/license-AGPL-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

## Introduction

Module for MSA management. Exportable functions in this module start with the prefix `pllmod_msa_`.

## Compilation instructions

Read the compilation instructions at the main README file

## Code

The code is written in C.

|    File              | Description                   |
|----------------------|-------------------------------|
|**pll_msa.c**         | Functions for analyzing MSAs. |

## Type definitions

* struct `pllmod_msa_stats_t`

## Functions

* `double * pllmod_msa_empirical_frequencies`
* `double * pllmod_msa_empirical_subst_rates`
* `double pllmod_msa_empirical_invariant_sites`
* `pllmod_msa_stats_t * pllmod_msa_compute_stats`
* `void pllmod_msa_destroy_stats`
* `pll_msa_t * pllmod_msa_filter`
* `pll_msa_t ** pllmod_msa_split`
* `int pllmod_msa_save_phylip`
