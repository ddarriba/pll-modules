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

struct default_params {
  pll_partition_t * partition;
  pll_unode_t * tree;
  unsigned int * params_indices;     /* indices for computing the likelihood */
  int gamma_mode;    /* discrete GAMMA rates computation mode (mean, median) */
};

struct freqs_params {
  pll_partition_t * partition;
  pll_unode_t * tree;
  unsigned int * params_indices;    /* indices for computing the likelihood */
  unsigned int params_index;        /* index of the frequencies to optimize */
  unsigned int highest_freq_state;  /* index of the highest frequency */
};

struct algo_subst_params {
  pll_partition_t * partition;
  pll_unode_t * tree;
  int * symmetries;                 /* substitution parameters symmetries */
  int subst_free_params;            /* number of free parameters */
  unsigned int * params_indices;    /* indices for computing the likelihood */
  unsigned int params_index;        /* index of the parameter to optimize */
};

struct rate_weights_params {
  pll_partition_t * partition;
  pll_unode_t * tree;
  unsigned int * params_indices;     /* indices for computing the likelihood */
  unsigned int highest_weight_state; /* index of the highest weight */
};

struct brlen_scaler_params {
  pll_partition_t * partition;
  pll_unode_t * tree;
  unsigned int * params_indices;     /* indices for computing the likelihood */
  double old_scaler;                 /* previous value of branch length scaler*/
};

struct treeinfo_opt_params {
  pllmod_treeinfo_t * treeinfo;
  int param_to_optimize;            /* which parameter is being optimized */
  unsigned int num_opt_partitions;  /* number of partitions to optimize */
  unsigned int params_index;        /* which matrix to optimize */
  unsigned int * num_free_params;   /* number of free params for each partition*/
  unsigned int * fixed_var_index;   /* which variable is not being optimized */
  treeinfo_param_set_cb param_set_cb;
};


/* optimize frequencies */
double target_freqs_func(void *p, double *x);

/* optimize substitution rates  */
double target_subst_params_func(void *p, double *x);

/* optimize alpha and discrete rate categories */
double target_alpha_func(void *p, double x);

/* optimize proportion of invariant sites */
double target_pinv_func(void *p, double x);

/* optimize alpha & proportion of invariant sites simultaneously */
double target_alpha_pinv_func(void *p, double *x);

/* optimize free rates */
double target_rates_func(void *p, double *x);

/* optimize free weights */
double target_weights_func(void *p, double *x);

/* optimize branch length scaler */
double target_brlen_scaler_func(void *p, double x);


double target_func_onedim_treeinfo(void *p, double *x, double *fx,
                                   int * converged);

double target_func_multidim_treeinfo(void * p, double ** x, double * fx,
                                     int * converged);

double target_subst_params_func_multi(void * p, double ** x, double * fx,
                                      int * converged);

double target_freqs_func_multi(void * p, double ** x, double * fx,
                               int * converged);


#endif /* ALGO_CALLBACK_H_ */
