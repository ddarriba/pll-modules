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
                          unsigned int * fixed_weight_index,
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

  /* find first state with frequency > min_freq, and use it as fixed freq */
  opt_params.fixed_freq_state = states;
  for (i = 0; i < states; ++i)
  {
    if (frequencies[i] > PLLMOD_OPT_MIN_FREQ)
    {
      opt_params.fixed_freq_state = i;
      break;
    }
  }

  assert(opt_params.fixed_freq_state < states);
  assert(frequencies[opt_params.fixed_freq_state] > 0.);

  cur_index = 0;
  for (i = 0; i < states; ++i)
  {
    if (i != opt_params.fixed_freq_state)
    {
      x[cur_index] = frequencies[i]
        / frequencies[opt_params.fixed_freq_state];
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
  lb[0] = min_alpha > 0. ? min_alpha : PLLMOD_OPT_MIN_ALPHA;
  ub[0] = max_alpha > 0. ? max_alpha : PLLMOD_OPT_MAX_ALPHA;
  bt[0] = PLLMOD_OPT_LBFGSB_BOUND_BOTH;

  /* init p-inv */
  x[1] = partition->prop_invar[params_indices[0]];
  lb[1] = min_pinv > PLLMOD_ALGO_LBFGSB_ERROR ? min_pinv :
      PLLMOD_OPT_MIN_PINV + PLLMOD_ALGO_LBFGSB_ERROR;
  ub[1] = max_pinv > 0. ? max_pinv : PLLMOD_OPT_MAX_PINV;
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

    fill_weights(weights, &(opt_params.fixed_weight_state), x,
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
                          unsigned int * fixed_weight_index,
                          double *x,
                          int *bt,
                          double *lb,
                          double *ub,
                          unsigned int n_weights)
{
  unsigned int i, cur_index = 0;

  *fixed_weight_index = n_weights;
  for (i = 0; i < n_weights; ++i)
  {
    if (weights[i] > PLLMOD_OPT_MIN_FREQ)
    {
      *fixed_weight_index = i;
      break;
    }
  }

  assert(*fixed_weight_index < n_weights);
  assert(weights[*fixed_weight_index] > 0.);

  for (i = 0; i < n_weights; ++i)
  {
    if (i != *fixed_weight_index)
    {
      bt[cur_index] = PLLMOD_OPT_LBFGSB_BOUND_BOTH;

      double r = weights[i] / weights[*fixed_weight_index];
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

static void fix_brlen_scalers(pllmod_treeinfo_t * treeinfo,
                              double min_scaler,
                              double max_scaler)
{
  unsigned int i;

  assert(treeinfo->brlen_scalers);

  double lowest_scaler  = treeinfo->brlen_scalers[0];
  double highest_scaler = treeinfo->brlen_scalers[0];

  /* collect brlen scaler for all partitions */
  if (treeinfo->parallel_reduce_cb)
  {
    for (i = 0; i < treeinfo->partition_count; ++i)
    {
      if (!treeinfo->partitions[i])
        treeinfo->brlen_scalers[i] = 0.0;
    }

    treeinfo->parallel_reduce_cb(treeinfo->parallel_context,
                                 treeinfo->brlen_scalers,
                                 treeinfo->partition_count,
                                 PLLMOD_COMMON_REDUCE_MAX);
  }

  /* skip remote partitions and those without rates/weight optimization */
  for (i = 0; i < treeinfo->partition_count; ++i)
  {
    if (treeinfo->brlen_scalers[i] < lowest_scaler)
      lowest_scaler = treeinfo->brlen_scalers[i];
    if (treeinfo->brlen_scalers[i] > highest_scaler)
      highest_scaler = treeinfo->brlen_scalers[i];
  }

  /* check if some scaler are out of bounds */
  if (lowest_scaler < min_scaler || highest_scaler > max_scaler)
  {
    double global_scaler;

    assert(lowest_scaler >= min_scaler || highest_scaler <= max_scaler);

    /* compute correction factor */
    if (lowest_scaler < min_scaler)
      global_scaler = min_scaler / lowest_scaler;
    else if (highest_scaler > max_scaler)
      global_scaler = max_scaler / highest_scaler;
    else
      assert(0);

    /* force scalers into bounds by multiplying with correction factor */
    for (i = 0; i < treeinfo->partition_count; ++i)
      treeinfo->brlen_scalers[i] *= global_scaler;

    /* multiply all branches by the inverse to preserve likelihood */
    pllmod_treeinfo_scale_branches_all(treeinfo, 1.0 / global_scaler);
  }
}

static pll_bool_t fix_brlen_minmax(pllmod_treeinfo_t * treeinfo, double blmin, double blmax)
{
  pll_bool_t brlen_fixed = PLL_FALSE;
  for (unsigned int i = 0; i < treeinfo->subnode_count; ++i)
  {
    pll_unode_t * snode = treeinfo->subnodes[i];
    if (snode->length < blmin)
    {
      pllmod_treeinfo_set_branch_length(treeinfo, snode, blmin);
      brlen_fixed = PLL_TRUE;
    }
    else if (snode->length > blmax)
    {
      pllmod_treeinfo_set_branch_length(treeinfo, snode, blmax);
      brlen_fixed = PLL_TRUE;
    }
  }

  return brlen_fixed;
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
  unsigned int param_count = 0;
  unsigned int i;

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
    unsigned int j = 0;
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
double pllmod_algo_opt_brlen_scalers_treeinfo(pllmod_treeinfo_t * treeinfo,
                                              double min_scaler,
                                              double max_scaler,
                                              double min_brlen,
                                              double max_brlen,
                                              double lh_epsilon)
{
  unsigned int i, j;
  double old_loglh, loglh;
  double *old_scalers = NULL;
  double *old_brlen   = NULL;

  if (treeinfo->brlen_linkage != PLLMOD_COMMON_BRLEN_SCALED)
  {
    pllmod_set_error(PLL_ERROR_PARAM_INVALID,
                     "Branch length scaler optimization works only in "
                     "scaled branch length mode.");
    return (double) PLL_FAILURE;
  }

  old_loglh = pllmod_treeinfo_compute_loglh(treeinfo, 0);

  /* save old brlen scalers in case we will have to revert optimization */
  old_scalers = (double *) calloc(treeinfo->init_partition_count,
                                  sizeof(double));
  old_brlen   = (double *) calloc(treeinfo->tree->edge_count,
                                  sizeof(double));
  for (i = 0, j = 0; i < treeinfo->partition_count; ++i)
  {
    if (treeinfo->partitions[i])
      old_scalers[j++] = treeinfo->brlen_scalers[i];
  }
  assert(j == treeinfo->init_partition_count);

  memcpy(old_brlen, treeinfo->branch_lengths[0],
         treeinfo->tree->edge_count * sizeof(double));

  /* make sure all brlen scalers are between min_scaler and max_scaler */
  fix_brlen_scalers(treeinfo, min_scaler, max_scaler);

  loglh = -1 * pllmod_algo_opt_onedim_treeinfo(treeinfo,
                                          PLLMOD_OPT_PARAM_BRANCH_LEN_SCALER,
                                          min_scaler,
                                          max_scaler,
                                          lh_epsilon);

  /* normalize scalers and scale the branches accordingly */
  pllmod_treeinfo_normalize_brlen_scalers(treeinfo);

  /* check that all branch lengths are within bounds after normalization,
   * and correct them as needed */
  pll_bool_t brlen_fixed = fix_brlen_minmax(treeinfo, min_brlen, max_brlen);

  if (brlen_fixed)
  {
    loglh = pllmod_treeinfo_compute_loglh(treeinfo, 0);
    if (loglh < old_loglh)
    {
      /* revert optimization and restore old values */
      for (i = 0, j = 0; i < treeinfo->partition_count; ++i)
      {
        if (treeinfo->partitions[i])
          treeinfo->brlen_scalers[i] = old_scalers[j++];
      }
      assert(j == treeinfo->init_partition_count);

      /* restore branch lengths */
      for (i = 0; i < treeinfo->subnode_count; ++i)
      {
        pll_unode_t * snode = treeinfo->subnodes[i];
        if (snode->node_index < snode->back->node_index)
        {
          pllmod_treeinfo_set_branch_length(treeinfo, snode,
                                            old_brlen[snode->pmatrix_index]);
        }
      }

      loglh = pllmod_treeinfo_compute_loglh(treeinfo, 0);
    }
  }

  free(old_scalers);
  free(old_brlen);

  return -1 * loglh;
}


PLL_EXPORT
double pllmod_algo_opt_subst_rates_treeinfo (pllmod_treeinfo_t * treeinfo,
                                             unsigned int params_index,
                                             double min_rate,
                                             double max_rate,
                                             double bfgs_factor,
                                             double tolerance)
{
  unsigned int i, j, k, l;

  const double factor = bfgs_factor > 0. ? bfgs_factor : PLLMOD_ALGO_BFGS_FACTR;

  double cur_logl;
  double **x, **lb, **ub;
  int **bt;

  unsigned int * subst_free_params;

  unsigned int part_count = 0;
  unsigned int max_free_params = 0;

  /* check how many alphas have to be optimized */
  for (i = 0; i < treeinfo->partition_count; ++i)
  {
    if (treeinfo->params_to_optimize[i] & PLLMOD_OPT_PARAM_SUBST_RATES)
    {
      part_count++;

      /* remote partition -> skip */
      if (!treeinfo->partitions[i])
        continue;

      unsigned int nrates = pllmod_util_subst_rate_count(treeinfo->partitions[i]->states);
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
    max_free_params = (unsigned int) tmp;
  }

  /* nothing to optimize */
  if (!part_count)
    return -1 * pllmod_treeinfo_compute_loglh(treeinfo, 0);

  x  = (double **) malloc(sizeof(double*) * (part_count));
  lb = (double **) malloc(sizeof(double*) * (part_count));
  ub = (double **) malloc(sizeof(double*) * (part_count));
  bt = (int **)    malloc(sizeof(int*)    * (part_count));
  subst_free_params = (unsigned int *) calloc(sizeof(unsigned int), part_count);

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

  unsigned int part = 0;
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

  unsigned int i, j;

  double cur_logl;
  double **x, **lb, **ub;
  int **bt;
  unsigned int * num_free_params;

  unsigned int part_count = 0;
  unsigned int max_free_params = 0;

  /* check how many frequencies have to be optimized */
  for (i = 0; i < treeinfo->partition_count; ++i)
  {
    if (treeinfo->params_to_optimize[i] & PLLMOD_OPT_PARAM_FREQUENCIES)
    {
      part_count++;

      /* remote partition -> skip */
      if (!treeinfo->partitions[i])
        continue;

      unsigned int nfree_params = treeinfo->partitions[i]->states - 1;
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
    max_free_params = (unsigned int) tmp;
  }

  x  = (double **) malloc(sizeof(double*) * part_count);
  lb = (double **) malloc(sizeof(double*) * part_count);
  ub = (double **) malloc(sizeof(double*) * part_count);
  bt = (int **)    malloc(sizeof(int*)    * part_count);
  num_free_params = (unsigned int *) calloc(sizeof(unsigned int), part_count);

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
  unsigned int part = 0;
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

#ifdef DEBUG
      printf("INITIAL freqs: ");
      for (size_t i = 0; i < states; ++i)
        printf("%f ", frequencies[i]);
      printf("\n");
#endif

    /* find first state with frequency > min_freq, and use it as fixed freq */
    unsigned int fixed_freq_state = states;
    for (j = 0; j < states; ++j)
    {
      if (frequencies[j] > min_freq)
      {
        fixed_freq_state = j;
        break;
      }
    }

    assert(fixed_freq_state < states);
    assert(frequencies[fixed_freq_state] > 0.);

    cur_index = 0;
    for (j = 0; j < states; ++j)
    {
      if (j != fixed_freq_state)
      {
        x[part][cur_index] = frequencies[j] / frequencies[fixed_freq_state];
        cur_index++;
      }
    }

    assert(cur_index == num_free_params[part]);

#ifdef DEBUG
      printf("INITIAL denorm freqs: ");
      for (size_t i = 0; i < cur_index; ++i)
        printf("%f ", x[part][i]);
      printf("\n");
#endif

    opt_params.fixed_var_index[part] = fixed_freq_state;

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

  unsigned int i;

  double cur_logl;
  double **x, **lb, **ub;
  double *xd;
  int **bt;
  unsigned int * num_free_params;

  unsigned int part_count = 0;
  unsigned int max_free_params = 2;

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
  num_free_params = (unsigned int *) calloc(sizeof(unsigned int), part_count);

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

  /* scale branch lengths such that likelihood is conserved */
  if (treeinfo->partition_count == 1)
    pllmod_treeinfo_scale_branches_all(treeinfo, brlen_scaler);
  else if (treeinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED)
    pllmod_treeinfo_scale_branches_partition(treeinfo, part_num, brlen_scaler);
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
  unsigned int i,j;

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

static void renormalize_free_rates(pllmod_treeinfo_t * treeinfo)
{
  unsigned int i, j;

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
}

PLL_EXPORT
double pllmod_algo_opt_rates_weights_treeinfo (pllmod_treeinfo_t * treeinfo,
                                               double min_rate,
                                               double max_rate,
                                               double min_brlen,
                                               double max_brlen,
                                               double bfgs_factor,
                                               double tolerance)
{
  const double factor = bfgs_factor > 0. ? bfgs_factor : PLLMOD_ALGO_BFGS_FACTR;

  unsigned int i;
  double old_logl, cur_logl, prev_logl;
  double **x, **lb, **ub;
  double *old_weights, *old_rates, *old_brlens, *old_scalers;
  int **bt;
  unsigned int * num_free_params;
  size_t rw_span;

  unsigned int part_count = 0;
  unsigned int local_part_count = 0;
  unsigned int part = 0;
  unsigned int local_part = 0;
  unsigned int max_free_params = 0;

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

      local_part_count++;

      unsigned int nfree_params = treeinfo->partitions[i]->rate_cats;
      if (nfree_params > max_free_params)
        max_free_params = nfree_params;
    }
  }

  old_logl = pllmod_treeinfo_compute_loglh(treeinfo, 0);

  /* nothing to optimize */
  if (!part_count)
    return -1 * old_logl;

  /* IMPORTANT: we need to know max_free_params among all threads! */
  if (treeinfo->parallel_reduce_cb)
  {
    double tmp = (double) max_free_params;
    treeinfo->parallel_reduce_cb(treeinfo->parallel_context, &tmp, 1,
                                 PLLMOD_COMMON_REDUCE_MAX);
    max_free_params = (unsigned int) tmp;
  }

  rw_span = max_free_params * sizeof(double);

  x  = (double **) calloc(sizeof(double*),  part_count);
  lb = (double **) calloc(sizeof(double*),  part_count);
  ub = (double **) calloc(sizeof(double*),  part_count);
  bt = (int **)    calloc(sizeof(int*),     part_count);
  num_free_params = (unsigned int *) calloc(sizeof(unsigned int), part_count);

  /* those values are the same for all partitions */
  lb[0] = (double *) malloc(sizeof(double) * (max_free_params));
  ub[0] = (double *) malloc(sizeof(double) * (max_free_params));
  bt[0] = (int *)    malloc(sizeof(int)    * (max_free_params));

  /* save old state for rollback: rates+weights+brlens+BL scalers */
  old_rates = old_weights = old_brlens = old_scalers = NULL;
  if (treeinfo->brlen_linkage != PLLMOD_COMMON_BRLEN_UNLINKED)
  {
    old_rates  = (double *) calloc(sizeof(double),  local_part_count * max_free_params);
    old_weights  = (double *) calloc(sizeof(double),  local_part_count * max_free_params);

    old_brlens  = (double *) calloc(sizeof(double),  treeinfo->tree->edge_count);
    memcpy(old_brlens, treeinfo->linked_branch_lengths, sizeof(double)*treeinfo->tree->edge_count);

    if (treeinfo->brlen_scalers)
    {
      old_scalers  = (double *) calloc(sizeof(double), treeinfo->partition_count);
      memcpy(old_scalers, treeinfo->brlen_scalers, sizeof(double)*treeinfo->partition_count);
    }
  }

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

      if (old_rates)
      {
        size_t rw_size = treeinfo->partitions[i]->rate_cats * sizeof(double);
        memcpy(old_rates + local_part * rw_span, treeinfo->partitions[i]->rates, rw_size);
        memcpy(old_weights + local_part * rw_span, treeinfo->partitions[i]->rate_weights, rw_size);
      }
      local_part++;
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
  renormalize_free_rates(treeinfo);

  /* update pmatrices and partials according to the new branches */
  cur_logl = pllmod_treeinfo_compute_loglh(treeinfo, 0);

  /* normalize scalers and scale the branches accordingly */
  if (treeinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_SCALED && treeinfo->partition_count > 1)
    pllmod_treeinfo_normalize_brlen_scalers(treeinfo);

  if (treeinfo->brlen_linkage != PLLMOD_COMMON_BRLEN_UNLINKED)
  {
    /* check that all branch lengths are within bounds after normalization,
     * and correct them as needed */
    pll_bool_t brlen_fixed = fix_brlen_minmax(treeinfo, min_brlen, max_brlen);

    if (brlen_fixed)
    {
      /* update pmatrices and partials according to the new branches */
      double new_logl = pllmod_treeinfo_compute_loglh(treeinfo, 0);

      DBG("pllmod_algo_opt_rates_weights_treeinfo: BRLEN_FIXED old/opt/fixed: %.12lf / %.12lf / %.12lf\n",
          old_logl, cur_logl, new_logl);

      if (new_logl < old_logl)
      {
        /* loglh worse than initial after enforcing min/max brlens -> rollback */
        assert(old_rates && old_weights && old_brlens);

        /* restore initial rates & weights */
        local_part = 0;
        for (i = 0; i < treeinfo->partition_count; ++i)
        {
          if (treeinfo->partitions[i] &&
              (treeinfo->params_to_optimize[i] &
              (PLLMOD_OPT_PARAM_FREE_RATES | PLLMOD_OPT_PARAM_RATE_WEIGHTS)))
          {
            size_t rw_size = treeinfo->partitions[i]->rate_cats * sizeof(double);
            memcpy(treeinfo->partitions[i]->rates, old_rates + local_part * rw_span, rw_size);
            memcpy(treeinfo->partitions[i]->rate_weights, old_weights + local_part * rw_span, rw_size);
            local_part++;
          }
        }
        /* restore initial branch lengths */
        for (i = 0; i < treeinfo->subnode_count; ++i)
          treeinfo->subnodes[i]->length = old_brlens[treeinfo->subnodes[i]->pmatrix_index];

        /* restore initial branch length scalers */
        if (treeinfo->brlen_scalers)
        {
          assert(old_scalers);
          memcpy(treeinfo->brlen_scalers, old_scalers, sizeof(double)*treeinfo->partition_count);
        }

        renormalize_free_rates(treeinfo);

        cur_logl = pllmod_treeinfo_compute_loglh(treeinfo, 0);

        DBG("pllmod_algo_opt_rates_weights_treeinfo: ROLLBACK, loglh = %.12lf\n", cur_logl);
      }
      else
        cur_logl = new_logl;
    }
  }

  /* cleanup */
  for (i = 0; i < part_count; ++i)
  {
    if (x[i])
      free(x[i]);
  }

  free(old_rates);
  free(old_weights);
  free(old_brlens);
  free(old_scalers);

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

PLL_EXPORT
double pllmod_algo_opt_brlen_treeinfo(pllmod_treeinfo_t * treeinfo,
                                      double min_brlen,
                                      double max_brlen,
                                      double lh_epsilon,
                                      int max_iters,
                                      int opt_method,
                                      int radius)
{
  return pllmod_opt_optimize_branch_lengths_local_multi(treeinfo->partitions,
                                                        treeinfo->partition_count,
                                                        treeinfo->root,
                                                        treeinfo->param_indices,
                                                        treeinfo->deriv_precomp,
                                                        treeinfo->branch_lengths,
                                                        treeinfo->brlen_scalers,
                                                        min_brlen,
                                                        max_brlen,
                                                        lh_epsilon,
                                                        max_iters,
                                                        radius,
                                                        1,    /* keep_update */
                                                        opt_method,
                                                        treeinfo->brlen_linkage,
                                                        treeinfo->parallel_context,
                                                        treeinfo->parallel_reduce_cb
                                                        );

}
