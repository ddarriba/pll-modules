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
#ifndef ALGO_CALLBACK_H_
#define ALGO_CALLBACK_H_

#include "pllmod_algorithm.h"

struct freqs_params {
  pll_partition_t * partition;
  pll_utree_t * tree;
  unsigned int * params_indices;    /* indices for computing the likelihood */
  unsigned int params_index;        /* index of the frequencies to optimize */
  unsigned int highest_freq_state;  /* index of the highest frequency */
};

struct algo_subst_params {
  pll_partition_t * partition;
  pll_utree_t * tree;
  int * symmetries;                 /* substitution parameters symmetries */
  int subst_free_params;            /* number of free parameters */
  unsigned int * params_indices;    /* indices for computing the likelihood */
  unsigned int params_index;        /* index of the parameter to optimize */
};

struct rate_weights_params {
  pll_partition_t * partition;
  pll_utree_t * tree;
  unsigned int * params_indices;     /* indices for computing the likelihood */
  unsigned int highest_weight_state; /* index of the highest weight */
};

/* optimize frequencies */
double target_freqs_func(void *p, double *x);

/* optimize substitution rates  */
double target_subst_params_func(void *p, double *x);

/* optimize free rates */
double target_rates_func(void *p, double *x);

/* optimize free weights */
double target_weights_func(void *p, double *x);

#endif /* ALGO_CALLBACK_H_ */
