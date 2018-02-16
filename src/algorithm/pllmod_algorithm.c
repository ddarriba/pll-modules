/*
 Copyright (C) 2016 Diego Darriba, Alexey Kozlov

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
  * @file pllmod_algorithm.c
  *
  * @brief High level algorithms
  *
  * This module contains high level algorithms that depend on several additional
  * modules.
  *
  * @author Diego Darriba
  * @author Alexey Kozlov
  */

#include "pllmod_algorithm.h"
#include "algo_callback.h"
#include "../pllmod_common.h"

static void fill_rates   (double *rates,
                          double *x,
                          int *bt, double *lb, double *ub,
                          double min_rate, double max_rate,
                          unsigned int n_rates);

static void fill_weights (double *weights,
                          unsigned int * highest_weight_index,
                          double *x, int *bt, double *lb, double *ub,
                          unsigned int n_weights);


PLL_EXPORT double pllmod_algo_opt_frequencies (pll_partition_t * partition,
                                               pll_unode_t * tree,
                                               unsigned int params_index,
                                               unsigned int * params_indices,
                                               double bfgs_factor,
                                               double tolerance)
{
  double cur_logl;
  double *x, *lb, *ub;
  int *bt;
  unsigned int i;
  struct freqs_params opt_params;
  opt_params.partition      = partition;
  opt_params.tree           = tree;
  opt_params.params_indices = params_indices;
  opt_params.params_index   = params_index;

  double * frequencies = partition->frequencies[params_index];
  unsigned int states  = partition->states;
  unsigned int cur_index;

  const double factor = bfgs_factor > 0. ? bfgs_factor : PLLMOD_ALGO_BFGS_FACTR;

  x  = (double *) malloc(sizeof(double) * (states - 1));
  lb = (double *) malloc(sizeof(double) * (states - 1));
  ub = (double *) malloc(sizeof(double) * (states - 1));
  bt = (int *)    malloc(sizeof(int)    * (states - 1));

  /* find highest frequency */
  opt_params.highest_freq_state = 0;
  for (i = 1; i < states; ++i)
    if (frequencies[i] > frequencies[opt_params.highest_freq_state])
      opt_params.highest_freq_state = i;

  cur_index = 0;
  for (i = 0; i < states; ++i)
  {
    if (i != opt_params.highest_freq_state)
    {
      x[cur_index] = frequencies[i]
        / frequencies[opt_params.highest_freq_state];
      lb[cur_index] = PLLMOD_OPT_MIN_FREQ;
      ub[cur_index] = PLLMOD_OPT_MAX_FREQ;
      bt[cur_index] = PLLMOD_OPT_LBFGSB_BOUND_BOTH;
      cur_index++;
    }
  }

  cur_logl = pllmod_opt_minimize_lbfgsb(x, lb, ub, bt,
                                 states-1,
                                 factor, tolerance,
                                 (void *) &opt_params,
                                 &target_freqs_func);

  /* update frequencies */
  target_freqs_func((void *)&opt_params, x);

  free(x);
  free(lb);
  free(ub);
  free(bt);

  return cur_logl;
}

PLL_EXPORT double pllmod_algo_opt_subst_rates (pll_partition_t * partition,
                                               pll_unode_t * tree,
                                               unsigned int params_index,
                                               unsigned int * params_indices,
                                               int * symmetries,
                                               double min_rate,
                                               double max_rate,
                                               double bfgs_factor,
                                               double tolerance)
{
  double cur_logl;
  double *x, *lb, *ub;
  int *bt;
  unsigned int i, j, k;

  double *subst_rates    = partition->subst_params[params_index];
  unsigned int states    = partition->states;
  unsigned int subst_params = (states * (states-1)) / 2;
  unsigned int subst_free_params;

  const double factor = bfgs_factor > 0. ? bfgs_factor : PLLMOD_ALGO_BFGS_FACTR;

  if (!symmetries)
  {
    subst_free_params = subst_params - 1;
  }
  else
  {
    subst_free_params = 0;
    for (i=0; i<subst_params; ++i)
    {
     if ((unsigned int)symmetries[i] > subst_free_params)
     {
      /* check that symmetries vector is correctly formatted */
      assert((unsigned int)symmetries[i] == (subst_free_params+1));
      ++subst_free_params;
     }
    }
  }

  struct algo_subst_params opt_params;
  opt_params.partition         = partition;
  opt_params.tree              = tree;
  opt_params.params_index      = params_index;
  opt_params.params_indices    = params_indices;
  opt_params.symmetries        = symmetries;
  opt_params.subst_free_params = subst_free_params;

  x  = (double *) malloc(sizeof(double) * (subst_free_params));
  lb = (double *) malloc(sizeof(double) * (subst_free_params));
  ub = (double *) malloc(sizeof(double) * (subst_free_params));
  bt = (int *)    malloc(sizeof(int)    * (subst_free_params));

  k = 0;
  for (i = 0; i < subst_free_params; ++i)
  {
    bt[i] = PLLMOD_OPT_LBFGSB_BOUND_BOTH;
    lb[i] = min_rate;
    ub[i] = max_rate;

    if (symmetries)
    {
      if ((unsigned int)symmetries[subst_params-1] == k)
        ++k;

      for (j=0; j<subst_params; ++j)
      {
        if ((unsigned int)symmetries[j] == k)
        {
          x[i] = subst_rates[j];
          break;
        }
      }
      ++k;
    }
    else
    {
      x[i] = subst_rates[i];
    }

    if (!x[i])
    {
      /* initialize to interval center */
      x[i] = (min_rate + max_rate) / 2.0;
    }
    if (x[i] < min_rate)
    {
      /* set to lower bound */
      x[i] = min_rate;
    }
    else if (x[i] > max_rate)
    {
      /* set to upper bound */
      x[i] = max_rate;
    }
  }

  cur_logl = pllmod_opt_minimize_lbfgsb(x, lb, ub, bt,
                                 subst_free_params,
                                 factor, tolerance,
                                 (void *) &opt_params, target_subst_params_func);

  free(x);
  free(lb);
  free(ub);
  free(bt);

  return cur_logl;
}

PLL_EXPORT double pllmod_algo_opt_alpha (pll_partition_t * partition,
                                         pll_unode_t * tree,
                                         unsigned int * params_indices,
                                         double min_alpha,
                                         double max_alpha,
                                         double *alpha,
                                         double tolerance)
{
  double cur_logl;
  double f2x;
  double xres;

  struct default_params opt_params;
  opt_params.partition      = partition;
  opt_params.tree           = tree;
  opt_params.params_indices = params_indices;
  opt_params.gamma_mode     = PLL_GAMMA_RATES_MEAN;  // for now

  xres = pllmod_opt_minimize_brent(min_alpha, *alpha, max_alpha,
                            tolerance,
                            &cur_logl,
                            &f2x,
                            (void *) &opt_params,
                            &target_alpha_func);

  cur_logl = target_alpha_func(&opt_params, xres);
  *alpha = xres;

  return cur_logl;
}

PLL_EXPORT double pllmod_algo_opt_pinv (pll_partition_t * partition,
                                        pll_unode_t * tree,
                                        unsigned int * params_indices,
                                        double min_pinv,
                                        double max_pinv,
                                        double tolerance)
{
  double cur_logl;
  double f2x;
  double xres;
  double start_pinv;
  struct default_params opt_params;
  opt_params.partition      = partition;
  opt_params.tree           = tree;
  opt_params.params_indices = params_indices;
  start_pinv = partition->prop_invar[params_indices[0]];

  xres = pllmod_opt_minimize_brent(min_pinv,
                            start_pinv,
                            max_pinv,
                            tolerance,
                            &cur_logl,
                            &f2x,
                            (void *) &opt_params,
                            &target_pinv_func);

  cur_logl = target_pinv_func(&opt_params, xres);

  return cur_logl;
}

PLL_EXPORT double pllmod_algo_opt_alpha_pinv (pll_partition_t * partition,
                                              pll_unode_t * tree,
                                              unsigned int * params_indices,
                                              double min_alpha,
                                              double max_alpha,
                                              double *alpha,
                                              double min_pinv,
                                              double max_pinv,
                                              double bfgs_factor,
                                              double tolerance)
{
  double cur_logl;
  double x[2], lb[2], ub[2];
  int bt[2];

  const double factor = bfgs_factor > 0. ? bfgs_factor : PLLMOD_ALGO_BFGS_FACTR;

  struct default_params opt_params;
  opt_params.partition      = partition;
  opt_params.tree           = tree;
  opt_params.params_indices = params_indices;
  opt_params.gamma_mode     = PLL_GAMMA_RATES_MEAN;  // for now

  /* init alpha */
  x[0] = *alpha;
  lb[0] = min_alpha ? min_alpha : PLLMOD_OPT_MIN_ALPHA;
  ub[0] = max_alpha ? max_alpha : PLLMOD_OPT_MAX_ALPHA;
  bt[0] = PLLMOD_OPT_LBFGSB_BOUND_BOTH;

  /* init p-inv */
  x[1] = partition->prop_invar[params_indices[0]];
  lb[1] = min_pinv > PLLMOD_ALGO_LBFGSB_ERROR ? min_pinv :
      PLLMOD_OPT_MIN_PINV + PLLMOD_ALGO_LBFGSB_ERROR;
  ub[1] = max_pinv ? max_pinv : PLLMOD_OPT_MAX_PINV;
  bt[1] = PLLMOD_OPT_LBFGSB_BOUND_BOTH;

  cur_logl = pllmod_opt_minimize_lbfgsb(x, lb, ub, bt,
                                 2,
                                 factor, tolerance,
                                 (void *) &opt_params,
                                 &target_alpha_pinv_func);

  /* save optimal alpha (p-inv is stored in the partition) */
  *alpha = x[0];

  return cur_logl;
}

PLL_EXPORT double pllmod_algo_opt_brlen_scaler (pll_partition_t * partition,
                                                pll_unode_t * root,
                                                unsigned int * params_indices,
                                                double * scaler,
                                                double min_scaler,
                                                double max_scaler,
                                                double tolerance)
{
  double cur_logl;
  double f2x;
  double xres;
  struct brlen_scaler_params opt_params;

  /* create a temporary tree with the scaled branches */
  pll_unode_t * scaled_tree = pll_utree_graph_clone(root);
  pllmod_utree_scale_branches_all(scaled_tree, *scaler);

  opt_params.partition      = partition;
  opt_params.tree           = scaled_tree;
  opt_params.params_indices = params_indices;
  opt_params.old_scaler     = *scaler;

  xres = pllmod_opt_minimize_brent(min_scaler,
                                   *scaler,
                                   max_scaler,
                                   tolerance,
                                   &cur_logl,
                                   &f2x,
                                   (void *) &opt_params,
                                   &target_brlen_scaler_func);

  cur_logl = target_brlen_scaler_func(&opt_params, xres);

  pll_utree_graph_destroy(scaled_tree, NULL);

  *scaler = xres;

  return cur_logl;
}

PLL_EXPORT double pllmod_algo_opt_rates_weights (pll_partition_t * partition,
                                                 pll_unode_t * tree,
                                                 unsigned int * params_indices,
                                                 double min_rate,
                                                 double max_rate,
                                                 double bfgs_factor,
                                                 double tolerance,
                                                 double * brlen_scaler,
                                                 int scale_branches)
{
  double cur_logl, prev_logl;
  double sum_weightrates, rate_scaler;
  double *x, *lb, *ub;
  int *bt;
  unsigned int i;

  double *rates          = partition->rates;
  double *weights        = partition->rate_weights;
  unsigned int rate_cats = partition->rate_cats;

  struct rate_weights_params opt_params;
  opt_params.partition      = partition;
  opt_params.tree           = tree;
  opt_params.params_indices = params_indices;

  const double factor = bfgs_factor > 0. ? bfgs_factor : PLLMOD_ALGO_BFGS_FACTR;

  x  = (double *) malloc(sizeof(double) * (rate_cats));
  lb = (double *) malloc(sizeof(double) * (rate_cats));
  ub = (double *) malloc(sizeof(double) * (rate_cats));
  bt = (int *)    malloc(sizeof(int)    * (rate_cats));

  /* 2 step BFGS */

  cur_logl = 0;
  do
  {
    prev_logl = cur_logl;

    /* optimize mixture weights */

    fill_weights(weights, &(opt_params.highest_weight_state), x,
                      bt, lb, ub, rate_cats);

    cur_logl = 1
        * pllmod_opt_minimize_lbfgsb (x, lb, ub, bt, rate_cats-1,
                                      factor, tolerance,
                                      (void *) &opt_params,
                                      target_weights_func);

    /* optimize mixture rates */

    fill_rates (rates,
                x, bt, lb, ub,
                min_rate, max_rate,
                rate_cats);

    cur_logl = pllmod_opt_minimize_lbfgsb(x, lb, ub, bt,
                                          rate_cats,
                                          factor, tolerance,
                                          (void *) &opt_params,
                                          target_rates_func);

  } while (!prev_logl || prev_logl - cur_logl > tolerance);

  /* force constraint sum(weights x rates) = 1.0 */
  sum_weightrates = 0.0;
  for (i=0; i<rate_cats; ++i)
    sum_weightrates += rates[i] * weights[i];
  rate_scaler = 1.0 / sum_weightrates;

  for (i=0; i<rate_cats; ++i)
    rates[i] *= rate_scaler;

  *brlen_scaler = sum_weightrates;

  if (scale_branches)
  {
    /* scale branch lengths such that likelihood is conserved */
    pllmod_utree_scale_branches_all(tree, sum_weightrates);

    /* update pmatrices and partials according to the new branches */
    cur_logl = -1 *
               pllmod_utree_compute_lk(partition,
                                       tree,
                                       params_indices,
                                       1,   /* update pmatrices */
                                       1);  /* update partials */
  }

  free(x);
  free(lb);
  free(ub);
  free(bt);

  return cur_logl;
}

/* STATIC FUNCTIONS */

static void fill_rates   (double *rates,
                          double *x,
                          int *bt, double *lb, double *ub,
                          double min_rate,
                          double max_rate,
                          unsigned int n_rates)
{
  unsigned int i;

  assert (min_rate > 1e-4 && max_rate > min_rate);

  for (i = 0; i < n_rates; ++i)
  {
    bt[i] = PLLMOD_OPT_LBFGSB_BOUND_BOTH;
    lb[i] = min_rate;
    ub[i] = max_rate;

    if (rates[i] < min_rate)
      x[i] = min_rate;
    else if (rates[i] > max_rate)
      x[i] = max_rate;
    else
      x[i] = rates[i];
  }
}

static void fill_weights (double *weights,
                          unsigned int * highest_weight_index,
                          double *x,
                          int *bt,
                          double *lb,
                          double *ub,
                          unsigned int n_weights)
{
  unsigned int i, cur_index = 0;

  *highest_weight_index = 0;
  for (i = 1; i < n_weights; ++i)
    if (weights[i] > weights[*highest_weight_index])
      *highest_weight_index = i;

  for (i = 0; i < n_weights; ++i)
  {
    if (i != *highest_weight_index)
    {
      bt[cur_index] = PLLMOD_OPT_LBFGSB_BOUND_BOTH;

      double r = weights[i] / weights[*highest_weight_index];
      lb[cur_index] = PLLMOD_ALGO_MIN_WEIGHT_RATIO;
      ub[cur_index] = PLLMOD_ALGO_MAX_WEIGHT_RATIO;
      if (r < lb[cur_index])
        x[cur_index] = lb[cur_index];
      else if (r > ub[cur_index])
        x[cur_index] = ub[cur_index];
      else
        x[cur_index] = r;

      cur_index++;
    }
  }
}

static int treeinfo_get_alpha(const pllmod_treeinfo_t * treeinfo,
                              unsigned int  part_num,
                              double * param_vals,
                              unsigned int param_count)
{
  if (part_num >= treeinfo->partition_count)
    return PLL_FAILURE;

  param_vals[0] = treeinfo->alphas[part_num];
  return PLL_SUCCESS;
}

static int treeinfo_set_alpha(pllmod_treeinfo_t * treeinfo,
                              unsigned int  part_num,
                              const double * param_vals,
                              unsigned int param_count)
{
  if (part_num >= treeinfo->partition_count)
    return PLL_FAILURE;

  treeinfo->alphas[part_num] = param_vals[0];

  pll_partition_t * partition = treeinfo->partitions[part_num];

  /* update rate categories */
  if (!pll_compute_gamma_cats (treeinfo->alphas[part_num],
                               partition->rate_cats,
                               partition->rates,
                               treeinfo->gamma_mode[part_num]))
    return PLL_FAILURE;

  return PLL_SUCCESS;
}

static int treeinfo_get_pinv(const pllmod_treeinfo_t * treeinfo,
                             unsigned int  part_num,
                             double * param_vals,
                             unsigned int param_count)
{
  if (part_num >= treeinfo->partition_count)
    return PLL_FAILURE;

  pll_partition_t * partition = treeinfo->partitions[part_num];
  param_vals[0] = partition->prop_invar[treeinfo->param_indices[part_num][0]];
  return PLL_SUCCESS;
}

static int treeinfo_set_pinv(pllmod_treeinfo_t * treeinfo,
                             unsigned int  part_num,
                             const double * param_vals,
                             unsigned int param_count)
{
  if (part_num >= treeinfo->partition_count)
    return PLL_FAILURE;

  unsigned int k;
  pll_partition_t * partition = treeinfo->partitions[part_num];

  /* update proportion of invariant sites */
  for (k = 0; k < partition->rate_cats; ++k)
  {
    if (!pll_update_invariant_sites_proportion(partition,
                                          treeinfo->param_indices[part_num][k],
                                          param_vals[0]))
      return PLL_FAILURE;
  }

  return PLL_SUCCESS;
}

static int treeinfo_get_brlen_scaler(const pllmod_treeinfo_t * treeinfo,
                                     unsigned int part_num,
                                     double * param_vals,
                                     unsigned int param_count)
{
  if (part_num >= treeinfo->partition_count)
    return PLL_FAILURE;

  param_vals[0] = treeinfo->brlen_scalers[part_num];
  return PLL_SUCCESS;
}

static int treeinfo_set_brlen_scaler(pllmod_treeinfo_t * treeinfo,
                                     unsigned int part_num,
                                     const double * param_vals,
                                     unsigned int param_count)
{
  if (part_num >= treeinfo->partition_count)
    return PLL_FAILURE;

  treeinfo->brlen_scalers[part_num] = param_vals[0];

  return PLL_SUCCESS;
}

PLL_EXPORT
double pllmod_algo_opt_onedim_treeinfo_custom(pllmod_treeinfo_t * treeinfo,
                                              int param_to_optimize,
                                              treeinfo_param_get_cb params_getter,
                                              treeinfo_param_set_cb params_setter,
                                              double min_value,
                                              double max_value,
                                              double tolerance)
{
  size_t param_count = 0;
  size_t i;

  /* check how many partitions have to be optimized */
  for (i = 0; i < treeinfo->partition_count; ++i)
  {
    if (treeinfo->params_to_optimize[i] & param_to_optimize)
      param_count++;
  }

  if (param_count > 0)
  {
    double * param_vals = (double *) malloc(param_count * sizeof(double));
    int * opt_mask = (int *) calloc(param_count, sizeof(int));

    /* collect current values of parameters */
    size_t j = 0;
    for (i = 0; i < treeinfo->partition_count; ++i)
    {
      if (treeinfo->params_to_optimize[i] & param_to_optimize)
      {
        pll_partition_t * partition = treeinfo->partitions[i];

        /* remote partition -> skip */
        if (!partition)
        {
          j++;
          continue;
        }

        params_getter(treeinfo, i, &param_vals[j], 1);
        opt_mask[j] = 1;
        j++;
      }
    }
    assert(j == param_count);

    struct treeinfo_opt_params opt_params;
    opt_params.treeinfo           = treeinfo;
    opt_params.param_to_optimize  = param_to_optimize;
    opt_params.num_opt_partitions = param_count;
    opt_params.param_set_cb       = params_setter;

    /* run BRENT optimization for all partitions in parallel */
    int ret = pllmod_opt_minimize_brent_multi(param_count,
                                              opt_mask,
                                              &min_value, param_vals, &max_value,
                                              tolerance, param_vals,
                                              NULL, NULL,  /* fx, f2x */
                                              (void *) &opt_params,
                                              &target_func_onedim_treeinfo,
                                              1 /* global_range */
                                              );

    free(param_vals);
    free(opt_mask);

    if (ret != PLL_SUCCESS)
    {
      assert(pll_errno);
      return -INFINITY;
    }
  }

  double cur_logl = pllmod_treeinfo_compute_loglh(treeinfo, 0);

  return -1 * cur_logl;
}

PLL_EXPORT double pllmod_algo_opt_onedim_treeinfo(pllmod_treeinfo_t * treeinfo,
                                                  int param_to_optimize,
                                                  double min_value,
                                                  double max_value,
                                                  double tolerance)
{
  if (__builtin_popcount(param_to_optimize) > 1)
  {
    pllmod_set_error(PLL_ERROR_PARAM_INVALID,
                     "Multi-parameter optimization is not supported by the "
                     "pllmod_algo_opt_onedim_treeinfo() function!");
    return -INFINITY;
  }

  treeinfo_param_get_cb params_getter = NULL;
  treeinfo_param_set_cb params_setter = NULL;

  switch (param_to_optimize)
  {
    case PLLMOD_OPT_PARAM_ALPHA:
      params_getter = treeinfo_get_alpha;
      params_setter = treeinfo_set_alpha;
      break;
    case PLLMOD_OPT_PARAM_PINV:
      params_getter = treeinfo_get_pinv;
      params_setter = treeinfo_set_pinv;
      break;
    case PLLMOD_OPT_PARAM_BRANCH_LEN_SCALER:
      params_getter = treeinfo_get_brlen_scaler;
      params_setter = treeinfo_set_brlen_scaler;
      break;
    default:
      pllmod_set_error(PLL_ERROR_PARAM_INVALID, "Unsupported parameter: %d",
                       param_to_optimize);
      return -INFINITY;
  }

  assert(params_getter && params_setter);

  return pllmod_algo_opt_onedim_treeinfo_custom(treeinfo,
                                                param_to_optimize,
                                                params_getter,
                                                params_setter,
                                                min_value,
                                                max_value,
                                                tolerance);
}

PLL_EXPORT
double pllmod_algo_opt_subst_rates_treeinfo (pllmod_treeinfo_t * treeinfo,
                                             unsigned int params_index,
                                             double min_rate,
                                             double max_rate,
                                             double bfgs_factor,
                                             double tolerance)
{
  size_t i, j, k, l;

  const double factor = bfgs_factor > 0. ? bfgs_factor : PLLMOD_ALGO_BFGS_FACTR;

  double cur_logl;
  double **x, **lb, **ub;
  int **bt;

  unsigned int * subst_free_params;

  size_t part_count = 0;
  size_t max_free_params = 0;

  /* check how many alphas have to be optimized */
  for (i = 0; i < treeinfo->partition_count; ++i)
  {
    if (treeinfo->params_to_optimize[i] & PLLMOD_OPT_PARAM_SUBST_RATES)
    {
      part_count++;

      /* remote partition -> skip */
      if (!treeinfo->partitions[i])
        continue;

      size_t nrates = pllmod_util_subst_rate_count(treeinfo->partitions[i]->states);
      if (nrates > max_free_params)
        max_free_params = nrates;
    }
  }

  /* IMPORTANT: we need to know max_free_params among all threads! */
  if (treeinfo->parallel_reduce_cb)
  {
    double tmp = (double) max_free_params;
    treeinfo->parallel_reduce_cb(treeinfo->parallel_context, &tmp, 1,
                                 PLLMOD_COMMON_REDUCE_MAX);
    max_free_params = (size_t) tmp;
  }

  /* nothing to optimize */
  if (!part_count)
    return -1 * pllmod_treeinfo_compute_loglh(treeinfo, 0);

  x  = (double **) malloc(sizeof(double*) * (part_count));
  lb = (double **) malloc(sizeof(double*) * (part_count));
  ub = (double **) malloc(sizeof(double*) * (part_count));
  bt = (int **)    malloc(sizeof(int*)    * (part_count));
  subst_free_params = (unsigned int *) malloc(sizeof(unsigned int) * (part_count));

  /* those values are the same for all partitions */
  lb[0] = (double *) malloc(sizeof(double) * (max_free_params));
  ub[0] = (double *) malloc(sizeof(double) * (max_free_params));
  bt[0] = (int *)    malloc(sizeof(int)    * (max_free_params));

  for (k = 0; k < max_free_params; ++k)
  {
    bt[0][k] = PLLMOD_OPT_LBFGSB_BOUND_BOTH;
    lb[0][k] = min_rate;
    ub[0][k] = max_rate;
  }

  size_t part = 0;
  for (i = 0; i < treeinfo->partition_count; ++i)
  {
    /* skip partition where no rate optimization is needed */
    if (!(treeinfo->params_to_optimize[i] & PLLMOD_OPT_PARAM_SUBST_RATES))
      continue;

    /* remote partition -> skip */
    if (!treeinfo->partitions[i])
    {
      x[part] = NULL;
      part++;
      continue;
    }

    pll_partition_t * partition = treeinfo->partitions[i];
    double * subst_rates    = partition->subst_params[params_index];
    unsigned int states    = partition->states;
    unsigned int subst_params = pllmod_util_subst_rate_count(states);
    int * symmetries = treeinfo->subst_matrix_symmetries[i];

    if (!symmetries)
    {
      subst_free_params[part] = subst_params - 1;
    }
    else
    {
      subst_free_params[part] = 0;
      for (k=0; k<subst_params; ++k)
      {
        if ((unsigned int)symmetries[k] > subst_free_params[part])
        {
          /* check that symmetries vector is correctly formatted */
          assert((unsigned int)symmetries[k] == (subst_free_params[part]+1));
          ++subst_free_params[part];
        }
      }
    }

    x[part]  = (double *) malloc(sizeof(double) * (subst_free_params[part]));
    bt[part] = bt[0];
    lb[part] = lb[0];
    ub[part] = ub[0];

    l = 0;
    for (k = 0; k < subst_free_params[part]; ++k)
    {
      if (symmetries)
      {
        if ((unsigned int)symmetries[subst_params-1] == l)
          ++l;

        for (j=0; j<subst_params; ++j)
        {
          if ((unsigned int)symmetries[j] == l)
          {
            x[part][k] = subst_rates[j];
            break;
          }
        }
        ++l;
      }
      else
      {
        x[part][k] = subst_rates[k];
      }

      if (x[part][k] < min_rate)
        x[part][k] = min_rate;
      else if (x[part][k] > max_rate)
        x[part][k] = max_rate;
    }

    part++;
  }

  assert(part == part_count);

  struct treeinfo_opt_params opt_params;
  opt_params.treeinfo           = treeinfo;
  opt_params.num_opt_partitions = part_count;
  opt_params.params_index       = params_index;
  opt_params.num_free_params    = subst_free_params;
  opt_params.fixed_var_index    = NULL;

  cur_logl = pllmod_opt_minimize_lbfgsb_multi(part_count, x, lb, ub, bt,
                                              subst_free_params,
                                              max_free_params,
                                              factor, tolerance,
                                              (void *) &opt_params,
                                              target_subst_params_func_multi);

  /* cleanup */
  for (i = 0; i < part_count; ++i)
  {
    if(x[i])
      free(x[i]);
  }

  free(lb[0]);
  free(ub[0]);
  free(bt[0]);

  free(x);
  free(lb);
  free(ub);
  free(bt);
  free(subst_free_params);

  return cur_logl;
}

PLL_EXPORT
double pllmod_algo_opt_frequencies_treeinfo (pllmod_treeinfo_t * treeinfo,
                                             unsigned int params_index,
                                             double min_freq,
                                             double max_freq,
                                             double bfgs_factor,
                                             double tolerance)
{
  const double factor = bfgs_factor > 0. ? bfgs_factor : PLLMOD_ALGO_BFGS_FACTR;

  size_t i, j;

  double cur_logl;
  double **x, **lb, **ub;
  int **bt;
  unsigned int * num_free_params;

  size_t part_count = 0;
  size_t max_free_params = 0;

  /* check how many frequencies have to be optimized */
  for (i = 0; i < treeinfo->partition_count; ++i)
  {
    if (treeinfo->params_to_optimize[i] & PLLMOD_OPT_PARAM_FREQUENCIES)
    {
      part_count++;

      /* remote partition -> skip */
      if (!treeinfo->partitions[i])
        continue;

      size_t nfree_params = treeinfo->partitions[i]->states - 1;
      if (nfree_params > max_free_params)
        max_free_params = nfree_params;
    }
  }

  /* nothing to optimize */
  if (!part_count)
    return -1 * pllmod_treeinfo_compute_loglh(treeinfo, 0);

  /* IMPORTANT: we need to know max_free_params among all threads! */
  if (treeinfo->parallel_reduce_cb)
  {
    double tmp = (double) max_free_params;
    treeinfo->parallel_reduce_cb(treeinfo->parallel_context, &tmp, 1,
                                 PLLMOD_COMMON_REDUCE_MAX);
    max_free_params = (size_t) tmp;
  }

  x  = (double **) malloc(sizeof(double*) * part_count);
  lb = (double **) malloc(sizeof(double*) * part_count);
  ub = (double **) malloc(sizeof(double*) * part_count);
  bt = (int **)    malloc(sizeof(int*)    * part_count);
  num_free_params = (unsigned int *) malloc(sizeof(unsigned int) * (part_count));

  /* those values are the same for all partitions */
  lb[0] = (double *) malloc(sizeof(double) * (max_free_params));
  ub[0] = (double *) malloc(sizeof(double) * (max_free_params));
  bt[0] = (int *)    malloc(sizeof(int)    * (max_free_params));

  for (j = 0; j < max_free_params; ++j)
  {
    bt[0][j] = PLLMOD_OPT_LBFGSB_BOUND_BOTH;
    lb[0][j] = min_freq;
    ub[0][j] = max_freq;
  }

  struct treeinfo_opt_params opt_params;
  opt_params.treeinfo       = treeinfo;
  opt_params.num_opt_partitions = part_count;
  opt_params.params_index   = params_index;
  opt_params.fixed_var_index = (unsigned int *) calloc(part_count,
                                                       sizeof(unsigned int));

  /* now iterate over partitions and collect current frequencies values */
  size_t part = 0;
  for (i = 0; i < treeinfo->partition_count; ++i)
  {
    // skip partition where no freqs optimization is needed
    if (!(treeinfo->params_to_optimize[i] & PLLMOD_OPT_PARAM_FREQUENCIES))
      continue;

    /* skip remote partitions (will be handled by other threads) */
    if (!treeinfo->partitions[i])
    {
      x[part] = NULL;
      part++;
      continue;
    }

    pll_partition_t * partition = treeinfo->partitions[i];
    double * frequencies = partition->frequencies[params_index];
    unsigned int states    = partition->states;
    unsigned int cur_index;

    num_free_params[part] = states - 1;

    x[part]  = (double *) malloc(sizeof(double) * (num_free_params[part]));
    bt[part] = bt[0];
    lb[part] = lb[0];
    ub[part] = ub[0];

    /* find highest frequency */
    unsigned int highest_freq_state = 0;
    for (j = 1; j < states; ++j)
      if (frequencies[j] > frequencies[highest_freq_state])
        highest_freq_state = j;

    cur_index = 0;
    for (j = 0; j < states; ++j)
    {
      if (j != highest_freq_state)
      {
        x[part][cur_index] = frequencies[j] / frequencies[highest_freq_state];
        cur_index++;
      }
    }

    opt_params.fixed_var_index[part] = highest_freq_state;

    part++;
  }

  assert(part == part_count);

  cur_logl = pllmod_opt_minimize_lbfgsb_multi(part_count, x, lb, ub, bt,
                                              num_free_params,
                                              max_free_params,
                                              factor, tolerance,
                                              (void *) &opt_params,
                                              target_freqs_func_multi);

  /* cleanup */
  for (i = 0; i < part_count; ++i)
    free(x[i]);

  free(lb[0]);
  free(ub[0]);
  free(bt[0]);

  free(x);
  free(lb);
  free(ub);
  free(bt);
  free(num_free_params);

  free(opt_params.fixed_var_index);

  return cur_logl;
}

PLL_EXPORT
double pllmod_algo_opt_alpha_pinv_treeinfo(pllmod_treeinfo_t * treeinfo,
                                           unsigned int params_index,
                                           double min_alpha,
                                           double max_alpha,
                                           double min_pinv,
                                           double max_pinv,
                                           double bfgs_factor,
                                           double tolerance)
{
  const double factor = bfgs_factor > 0. ? bfgs_factor : PLLMOD_ALGO_BFGS_FACTR;
  const int params_to_optimize = PLLMOD_OPT_PARAM_ALPHA | PLLMOD_OPT_PARAM_PINV;

  size_t i;

  double cur_logl;
  double **x, **lb, **ub;
  double *xd;
  int **bt;
  unsigned int * num_free_params;

  size_t part_count = 0;
  size_t max_free_params = 2;

  /* check in how many partitions both alpha AND p-inv have to be optimized */
  for (i = 0; i < treeinfo->partition_count; ++i)
  {
    if ((treeinfo->params_to_optimize[i] & params_to_optimize) == params_to_optimize)
    {
      part_count++;
    }
  }

  /* nothing to optimize */
  if (!part_count)
    return -1 * pllmod_treeinfo_compute_loglh(treeinfo, 0);

  x  = (double **) malloc(sizeof(double*) * part_count);
  xd = (double *)  malloc(sizeof(double)  * part_count * 2);
  lb = (double **) malloc(sizeof(double*) * part_count);
  ub = (double **) malloc(sizeof(double*) * part_count);
  bt = (int **)    malloc(sizeof(int*)    * part_count);
  num_free_params = (unsigned int *) malloc(sizeof(unsigned int) * (part_count));

  /* those values are the same for all partitions */
  lb[0] = (double *) malloc(sizeof(double) * 2);
  ub[0] = (double *) malloc(sizeof(double) * 2);
  bt[0] = (int *)    malloc(sizeof(int)    * 2);

  /* init bounds for alpha & p-inv */
  lb[0][0] = min_alpha ? min_alpha : PLLMOD_OPT_MIN_ALPHA;
  ub[0][0] = max_alpha ? max_alpha : PLLMOD_OPT_MAX_ALPHA;
  bt[0][0] = PLLMOD_OPT_LBFGSB_BOUND_BOTH;

  lb[0][1] = min_pinv > PLLMOD_ALGO_LBFGSB_ERROR ? min_pinv :
      PLLMOD_OPT_MIN_PINV + PLLMOD_ALGO_LBFGSB_ERROR;
  ub[0][1] = max_pinv ? max_pinv : PLLMOD_OPT_MAX_PINV;
  bt[0][1] = PLLMOD_OPT_LBFGSB_BOUND_BOTH;

  struct treeinfo_opt_params opt_params;
  opt_params.treeinfo       = treeinfo;
  opt_params.param_to_optimize = params_to_optimize;
  opt_params.num_opt_partitions = part_count;
  opt_params.params_index   = params_index;
  opt_params.fixed_var_index = NULL;

  /* now iterate over partitions and collect current alpha/p-inv values */
  size_t part = 0;
  for (i = 0; i < treeinfo->partition_count; ++i)
  {
    // skip partition where no freqs optimization is needed
    if ((treeinfo->params_to_optimize[i] & params_to_optimize) != params_to_optimize)
      continue;

    /* skip remote partitions (will be handled by other threads) */
    if (!treeinfo->partitions[i])
    {
      x[part] = NULL;
      part++;
      continue;
    }

    pll_partition_t * partition = treeinfo->partitions[i];

    /* init alpha & p-inv */
    x[part]  = xd + part * 2;
    x[part][0] = treeinfo->alphas[i];
    x[part][1] = partition->prop_invar[params_index];
    num_free_params[part] = max_free_params;

    bt[part] = bt[0];
    lb[part] = lb[0];
    ub[part] = ub[0];

    part++;
  }

  assert(part == part_count);

  cur_logl = pllmod_opt_minimize_lbfgsb_multi(part_count, x, lb, ub, bt,
                                              num_free_params,
                                              max_free_params,
                                              factor, tolerance,
                                              (void *) &opt_params,
                                              target_func_multidim_treeinfo);

  /* cleanup */
  free(lb[0]);
  free(ub[0]);
  free(bt[0]);

  free(x);
  free(xd);
  free(lb);
  free(ub);
  free(bt);
  free(num_free_params);

  return cur_logl;
}

static void scales_rates_and_branches(pllmod_treeinfo_t * treeinfo,
                                      size_t part_num,
                                      double rate_scaler)
{
  assert(treeinfo);
  assert(part_num < treeinfo->partition_count);
  assert(rate_scaler > 0.);

  pll_partition_t * partition = treeinfo->partitions[part_num];
  double * rates              = partition->rates;
  unsigned int rate_cats      = partition->rate_cats;
  double brlen_scaler 		  = 1.0 / rate_scaler;
  size_t j;

  for (j = 0; j < rate_cats; ++j)
    rates[j] *= rate_scaler;

  const int scale_branches = (treeinfo->partition_count == 1 ||
		  	  	  	   treeinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED);

  if (scale_branches)
  {
    //TODO adapt for unlinked branches

    /* scale branch lengths such that likelihood is conserved */
    pllmod_utree_scale_branches_all(treeinfo->root, brlen_scaler);
  }
  else
  {
    assert(treeinfo->brlen_scalers);

    /* update brlen scalers */
      treeinfo->brlen_scalers[part_num] *= brlen_scaler;
  }
}

static void fix_free_rates(pllmod_treeinfo_t * treeinfo,
                           double min_rate,
                           double max_rate)
{
  size_t i,j;

  for (i = 0; i < treeinfo->partition_count; ++i)
  {
    /* skip remote partitions and those without rates/weight optimization */
    if (!treeinfo->partitions[i] ||
    	!(treeinfo->params_to_optimize[i] & PLLMOD_OPT_PARAM_FREE_RATES))
      continue;

    pll_partition_t * partition = treeinfo->partitions[i];
    double * rates              = partition->rates;
    unsigned int rate_cats      = partition->rate_cats;
    double lowest_rate 			= rates[0];
    double highest_rate 		= rates[0];
    double rate_scaler;

    /* force constraint sum(weights x rates) = 1.0 */
    for (j = 1; j < rate_cats; ++j)
    {
      if (rates[j] < lowest_rate)
        lowest_rate = rates[j];
      if (rates[j] > highest_rate)
        highest_rate = rates[j];
    }

    if (lowest_rate < min_rate || highest_rate > max_rate)
    {
      assert(lowest_rate >= min_rate || highest_rate <= max_rate);

      if (lowest_rate < min_rate)
        rate_scaler = min_rate / lowest_rate;
      else if (highest_rate > max_rate)
        rate_scaler = max_rate / highest_rate;
      else
        assert(0);

      scales_rates_and_branches(treeinfo, i, rate_scaler);
    }
  }
}

PLL_EXPORT
double pllmod_algo_opt_rates_weights_treeinfo (pllmod_treeinfo_t * treeinfo,
                                               double min_rate,
                                               double max_rate,
                                               double bfgs_factor,
                                               double tolerance)
{
  const double factor = bfgs_factor > 0. ? bfgs_factor : PLLMOD_ALGO_BFGS_FACTR;

  size_t i, j;
  double cur_logl, prev_logl;
  double **x, **lb, **ub;
  int **bt;
  unsigned int * num_free_params;

  size_t part_count = 0;
  size_t max_free_params = 0;

  /* check how many frequencies have to be optimized */
  for (i = 0; i < treeinfo->partition_count; ++i)
  {
    if (treeinfo->params_to_optimize[i] &
        (PLLMOD_OPT_PARAM_FREE_RATES | PLLMOD_OPT_PARAM_RATE_WEIGHTS))
    {
      part_count++;

      /* remote partition -> skip */
      if (!treeinfo->partitions[i])
        continue;

      size_t nfree_params = treeinfo->partitions[i]->rate_cats;
      if (nfree_params > max_free_params)
        max_free_params = nfree_params;
    }
  }

  /* nothing to optimize */
  if (!part_count)
    return -1 * pllmod_treeinfo_compute_loglh(treeinfo, 0);

  /* IMPORTANT: we need to know max_free_params among all threads! */
  if (treeinfo->parallel_reduce_cb)
  {
    double tmp = (double) max_free_params;
    treeinfo->parallel_reduce_cb(treeinfo->parallel_context, &tmp, 1,
                                 PLLMOD_COMMON_REDUCE_MAX);
    max_free_params = (size_t) tmp;
  }

  x  = (double **) calloc(sizeof(double*),  part_count);
  lb = (double **) calloc(sizeof(double*),  part_count);
  ub = (double **) calloc(sizeof(double*),  part_count);
  bt = (int **)    calloc(sizeof(int*),     part_count);
  num_free_params = (unsigned int *) malloc(sizeof(unsigned int) * (part_count));

  /* those values are the same for all partitions */
  lb[0] = (double *) malloc(sizeof(double) * (max_free_params));
  ub[0] = (double *) malloc(sizeof(double) * (max_free_params));
  bt[0] = (int *)    malloc(sizeof(int)    * (max_free_params));

  size_t part = 0;
  for (i = 0; i < treeinfo->partition_count; ++i)
  {
    if (!(treeinfo->params_to_optimize[i] &
        (PLLMOD_OPT_PARAM_FREE_RATES | PLLMOD_OPT_PARAM_RATE_WEIGHTS)))
      continue;

    if (treeinfo->partitions[i])
    {
      x[part]  = (double *) malloc(sizeof(double) * (max_free_params));
      lb[part] = lb[0];
      ub[part] = ub[0];
      bt[part] = bt[0];
    }
    part++;
  }
  assert(part == part_count);

  struct treeinfo_opt_params opt_params;
  opt_params.treeinfo       = treeinfo;
  opt_params.num_opt_partitions = part_count;
  opt_params.params_index   = 0;
  opt_params.fixed_var_index = (unsigned int *) calloc(part_count,
                                                       sizeof(unsigned int));

  /* check if we have rates which are outside the bounds, and correct them by scaling */
  fix_free_rates(treeinfo, min_rate, max_rate);

  /* 2 step BFGS */
  cur_logl = -1 * pllmod_treeinfo_compute_loglh(treeinfo, 0);
  DBG("pllmod_algo_opt_rates_weights_treeinfo: START: logLH = %.15lf\n", cur_logl);
  do
  {
    prev_logl = cur_logl;

    /* optimize mixture weights */
    size_t part = 0;
    for (i = 0; i < treeinfo->partition_count; ++i)
    {
      if (treeinfo->params_to_optimize[i] & PLLMOD_OPT_PARAM_RATE_WEIGHTS)
      {
        pll_partition_t * partition = treeinfo->partitions[i];

        /* remote partition -> skip */
        if (!partition)
        {
          part++;
          continue;
        }

        num_free_params[part] = partition->rate_cats - 1;

        fill_weights(partition->rate_weights,
                     &(opt_params.fixed_var_index[part]), x[part],
                     bt[part], lb[part], ub[part], partition->rate_cats);
        part++;
      }
    }

    assert(part == part_count);

    opt_params.param_to_optimize = PLLMOD_OPT_PARAM_RATE_WEIGHTS;

    cur_logl = pllmod_opt_minimize_lbfgsb_multi(part_count, x, lb, ub, bt,
                                                num_free_params,
                                                max_free_params,
                                                factor, tolerance,
                                                (void *) &opt_params,
                                                target_func_multidim_treeinfo);

    DBG("pllmod_algo_opt_rates_weights_treeinfo: AFTER WEIGHTS: logLH = %.15lf\n", cur_logl);

    /* optimize mixture rates */

    part = 0;
    for (i = 0; i < treeinfo->partition_count; ++i)
    {
      if (treeinfo->params_to_optimize[i] & PLLMOD_OPT_PARAM_FREE_RATES)
      {
        pll_partition_t * partition = treeinfo->partitions[i];

        /* remote partition -> skip */
        if (!partition)
        {
          part++;
          continue;
        }

        num_free_params[part] = partition->rate_cats;

        DBG("pllmod_algo_opt_rates_weights_treeinfo: OLD RATES = (%.12lf %.12lf %.12lf %.12lf)\n",
            partition->rates[0], partition->rates[1], partition->rates[2], partition->rates[3]);

        fill_rates (partition->rates,
                    x[part], bt[part], lb[part], ub[part],
                    min_rate, max_rate,
                    partition->rate_cats);

        part++;
      }
    }

    opt_params.param_to_optimize = PLLMOD_OPT_PARAM_FREE_RATES;

    cur_logl = pllmod_opt_minimize_lbfgsb_multi(part_count, x, lb, ub, bt,
                                                num_free_params,
                                                max_free_params,
                                                factor, tolerance,
                                                (void *) &opt_params,
                                                target_func_multidim_treeinfo);

    DBG("pllmod_algo_opt_rates_weights_treeinfo: AFTER RATES: logLH = %.15lf\n", cur_logl);
  }
  while (prev_logl - cur_logl > tolerance);

  /* now re-normalize rates and scale the branches accordingly */
  for (i = 0; i < treeinfo->partition_count; ++i)
  {
    /* skip remote partitions and those without rates/weight optimization */
    if (!treeinfo->partitions[i] || !(treeinfo->params_to_optimize[i] &
        (PLLMOD_OPT_PARAM_FREE_RATES | PLLMOD_OPT_PARAM_RATE_WEIGHTS)))
      continue;

    pll_partition_t * partition = treeinfo->partitions[i];
    double * rates              = partition->rates;
    double * weights            = partition->rate_weights;
    unsigned int rate_cats      = partition->rate_cats;
    double sum_weightrates, rate_scaler;

    /* force constraint sum(weights x rates) = 1.0 */
    sum_weightrates = 0.0;
    for (j = 0; j < rate_cats; ++j)
      sum_weightrates += rates[j] * weights[j];
    rate_scaler = 1.0 / sum_weightrates;


    scales_rates_and_branches(treeinfo, i, rate_scaler);
  }

  /* update pmatrices and partials according to the new branches */
  cur_logl = pllmod_treeinfo_compute_loglh(treeinfo, 0);

  /* cleanup */
  for (i = 0; i < part_count; ++i)
  {
    if (x[i])
      free(x[i]);
  }

  free(lb[0]);
  free(ub[0]);
  free(bt[0]);

  free(x);
  free(lb);
  free(ub);
  free(bt);
  free(num_free_params);

  free(opt_params.fixed_var_index);

  return -1 * cur_logl;
}
