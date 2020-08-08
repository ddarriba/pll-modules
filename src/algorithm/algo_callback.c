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
  * @file algo_callback.c
  *
  * @brief Callback functions for algorithms module
  *
  * @author Diego Darriba
  * @author Alexey Kozlov
  */

#include "../pllmod_common.h"
#include "algo_callback.h"

double target_freqs_func(void *p, double *x)
{
  struct freqs_params * params = (struct freqs_params *) p;

  pll_partition_t * partition     = params->partition;
  pll_unode_t * root              = params->tree;
  unsigned int * params_indices   = params->params_indices;
  unsigned int params_index       = params->params_index;
  unsigned int highest_freq_state = params->highest_freq_state;
  unsigned int states             = partition->states;
  double *freqs                   = partition->frequencies[params_index];
  double sum_ratios = 1.0;

  unsigned int i, cur_index;
  double score;

  /* update frequencies */
  for (i = 0; i < (states - 1); ++i)
  {
    assert(x[i] == x[i]);
    sum_ratios += x[i];
  }
  cur_index = 0;
  for (i = 0; i < states; ++i)
  {
    if (i != highest_freq_state)
    {
      freqs[i] = x[cur_index] / sum_ratios;
      cur_index++;
    }
  }
  freqs[highest_freq_state] = 1.0 / sum_ratios;

  /* important!! invalidate eigen-decomposition */
  partition->eigen_decomp_valid[params_index] = 0;

  /* compute negative score */
  score = -1 * pllmod_utree_compute_lk(partition,
                                       root,
                                       params_indices,
                                       1,   /* update pmatrices */
                                       1);  /* update partials */
  return score;
}

double target_subst_params_func(void *p, double *x)
{
  struct algo_subst_params * params = (struct algo_subst_params *) p;

  unsigned int i,j,k;
  pll_partition_t * partition     = params->partition;
  pll_unode_t * root              = params->tree;
  unsigned int * params_indices   = params->params_indices;
  unsigned int params_index       = params->params_index;
  unsigned int subst_free_params  = params->subst_free_params;
  int * symmetries                = params->symmetries;

  unsigned int states             = partition->states;
  unsigned int subst_params       = (states * (states-1))/2;
  double *subst_rates             = partition->subst_params[params_index];

  if (symmetries)
  {
    /* assign values to the substitution rates */
    k = 0;
    for (i = 0; i <= subst_free_params; ++i)
    {
      double next_value =
               (i == (unsigned int)symmetries[subst_params - 1]) ? 1.0 : x[k++];
      for (j = 0; j < subst_params; j++)
      {
        if ((unsigned int)symmetries[j] == i)
        {
          subst_rates[j] = next_value;
        }
      }
    }
  }
  else
  {
    memcpy (subst_rates, x, ((size_t)subst_params - 1) * sizeof(double));
  }

  /* important!! invalidate eigen-decomposition */
  partition->eigen_decomp_valid[params_index] = 0;

  /* compute negative score */
  double score = -1 *
                 pllmod_utree_compute_lk(partition,
                                         root,
                                         params_indices,
                                         1,   /* update pmatrices */
                                         1);  /* update partials */
  return score;
}

double target_alpha_func(void *p, double x)
{
  struct default_params * params = (struct default_params *) p;
  pll_partition_t * partition   = params->partition;
  pll_unode_t * root            = params->tree;
  unsigned int * params_indices = params->params_indices;

  /* update rate categories */
  if (!pll_compute_gamma_cats (x,
                               partition->rate_cats,
                               partition->rates,
                               params->gamma_mode))
  {
    return PLL_FAILURE;
  }

  /* compute negative score */
  double score = -1 *
                 pllmod_utree_compute_lk(partition,
                                         root,
                                         params_indices,
                                         1,   /* update pmatrices */
                                         1);  /* update partials */
  return score;
}

double target_pinv_func(void *p, double x)
{
  struct default_params * params = (struct default_params *) p;
  pll_partition_t * partition   = params->partition;
  pll_unode_t * root            = params->tree;
  unsigned int * params_indices = params->params_indices;
  unsigned int i;

  /* update proportion of invariant sites */
  for (i=0; i<partition->rate_cats; ++i)
    pll_update_invariant_sites_proportion(partition,
                                          params_indices[i],
                                          x);

  /* compute negative score */
  double score = -1 *
                 pllmod_utree_compute_lk(partition,
                                         root,
                                         params_indices,
                                         1,   /* update pmatrices */
                                         1);  /* update partials */
  return score;
}

double target_alpha_pinv_func(void *p, double *x)
{
  struct default_params * params = (struct default_params *) p;
  pll_partition_t * partition   = params->partition;
  pll_unode_t * root            = params->tree;
  unsigned int * params_indices = params->params_indices;
  unsigned int i;

  /* update rate categories */
  if (!pll_compute_gamma_cats (x[0],
                               partition->rate_cats,
                               partition->rates,
                               params->gamma_mode))
  {
    return PLL_FAILURE;
  }

  /* update proportion of invariant sites */
  for (i=0; i<partition->rate_cats; ++i)
    pll_update_invariant_sites_proportion(partition,
                                          params_indices[i],
                                          x[1]);

  /* compute negative score */
  double score = -1 *
                 pllmod_utree_compute_lk(partition,
                                         root,
                                         params_indices,
                                         1,   /* update pmatrices */
                                         1);  /* update partials */
  return score;
}

double target_rates_func(void *p, double *x)
{
  struct rate_weights_params * params = (struct rate_weights_params *) p;
  pll_partition_t * partition   = params->partition;
  pll_unode_t * root            = params->tree;
  unsigned int * params_indices = params->params_indices;

  /* update rate categories */
  memcpy(partition->rates, x, partition->rate_cats*sizeof(double));

  /* compute negative score */
  double score = -1 *
                 pllmod_utree_compute_lk(partition,
                                         root,
                                         params_indices,
                                         1,   /* update pmatrices */
                                         1);  /* update partials */
  return score;
}

double target_weights_func(void *p, double *x)
{
  struct rate_weights_params * params = (struct rate_weights_params *) p;

  pll_partition_t * partition       = params->partition;
  pll_unode_t * root                = params->tree;
  unsigned int * params_indices     = params->params_indices;
  unsigned int highest_weight_state = params->highest_weight_state;
  unsigned int n_weights            = partition->rate_cats;
  double sum_ratios = 1.0;

  double * weights = partition->rate_weights;
  unsigned int i, cur_weight;
  double score;

  for (i = 0; i < (n_weights - 1); ++i)
    sum_ratios += x[i];

  cur_weight = 0;
  for (i = 0; i < (n_weights); ++i)
    if (i != highest_weight_state)
    {
      weights[i] = x[cur_weight++] / sum_ratios;
    }
  weights[highest_weight_state] = 1.0 / sum_ratios;

  /* update weights */
  // memcpy(partition->rate_weights, weights, partition->rate_cats*sizeof(double));

  /* compute negative score */
  score = -1 * pllmod_utree_compute_lk(partition,
                                       root,
                                       params_indices,
                                       0,   /* update pmatrices */
                                       0);  /* update partials */

  return score;
}

double target_brlen_scaler_func(void *p, double x)
{
  struct brlen_scaler_params * params = (struct brlen_scaler_params *) p;
  pll_partition_t * partition   = params->partition;
  pll_unode_t * root            = params->tree;
  unsigned int * params_indices = params->params_indices;

  /* scale branches according to the new factor */
  pllmod_utree_scale_branches_all(root, x / params->old_scaler);

  /* store the old scaler value */
  params->old_scaler = x;

  /* compute negative score */
  double score = -1 *
                 pllmod_utree_compute_lk(partition,
                                         root,
                                         params_indices,
                                         1,   /* update pmatrices */
                                         1);  /* update partials */
  return score;
}

double target_func_onedim_treeinfo(void *p, double *x, double *fx, int * converged)
{
  struct treeinfo_opt_params * params = (struct treeinfo_opt_params *) p;

  pllmod_treeinfo_t * treeinfo        = params->treeinfo;
  int param_to_optimize               = params->param_to_optimize;
  unsigned int num_parts              = params->num_opt_partitions;
  treeinfo_param_set_cb param_setter  = params->param_set_cb;

  double score = -INFINITY;

  /* any partitions which have not converged yet? */
  double unconverged_flag = 0.;

  size_t i, j=0;
  for (i = 0; i < treeinfo->partition_count; ++i)
  {
    pll_partition_t * partition = treeinfo->partitions[i];

    if (treeinfo->params_to_optimize[i] & param_to_optimize)
    {
      if (!partition || (converged && converged[j]))
      {
        /* partitions has converged, skip it */
        j++;
        continue;
      }

      unconverged_flag = 1.;

      /* if x=NULL, function was called solely to check convergence
       * -> no update of parameter values & LH computation */
      if (x)
        param_setter(treeinfo, i, &x[j], 1);

      j++;
    }
  }

  assert(j == num_parts);

  /* compute negative score */
  if (x)
    score = -1 * pllmod_treeinfo_compute_loglh(treeinfo, 0);

//  printf("score: %lf\n", score);

  /* copy per-partition likelihood to the output array */
  if (fx)
  {
    j = 0;
    for (i = 0; i < treeinfo->partition_count; ++i)
      if (treeinfo->params_to_optimize[i] & param_to_optimize)
        fx[j++] = -1 * treeinfo->partition_loglh[i];
  }

  if (converged)
  {
    /* check if there is at least one unconverged partition in *any* thread */
    if (treeinfo->parallel_reduce_cb)
    {
      treeinfo->parallel_reduce_cb(treeinfo->parallel_context, &unconverged_flag, 1,
                                   PLLMOD_COMMON_REDUCE_SUM);
    }
    converged[num_parts] = unconverged_flag > 0. ? 0 : 1;
  }

  return score;
}

double target_func_multidim_treeinfo(void * p, double ** x, double * fx,
                                     int * converged)
{
  struct treeinfo_opt_params * params = (struct treeinfo_opt_params *) p;

  pllmod_treeinfo_t * treeinfo      = params->treeinfo;
  unsigned int num_parts            = params->num_opt_partitions;
  unsigned int * fixed_var_index    = params->fixed_var_index;
  int params_to_optimize            = params->param_to_optimize;

  double score = -INFINITY;

  /* any partitions which have not converged yet? */
  double unconverged_flag = 0.;

  size_t i, j;
  size_t part = 0;
  for (i = 0; i < treeinfo->partition_count; ++i)
  {
    pll_partition_t * partition = treeinfo->partitions[i];

    if ((treeinfo->params_to_optimize[i] & params_to_optimize) == params_to_optimize)
    {
      if (!partition || (converged && converged[part]))
      {
        /* partitions has converged, skip it */
        part++;
        continue;
      }

      unconverged_flag = 1.;

      /* function was called solely to check convergence -> no LH computation */
      if (!x)
      {
        part++;
        continue;
      }

      switch (params_to_optimize)
      {
        case PLLMOD_OPT_PARAM_ALPHA | PLLMOD_OPT_PARAM_PINV:
          /* update GAMMA rate categories */
          treeinfo->alphas[i] = x[part][0];
          if (!pll_compute_gamma_cats (treeinfo->alphas[i],
                                       partition->rate_cats,
                                       partition->rates,
                                       params->treeinfo->gamma_mode[i]))
          {
            assert(pll_errno);
            return PLL_FAILURE;
          }

          /* update proportion of invariant sites */
          for (j=0; j<partition->rate_cats; ++j)
          {
            if (!pll_update_invariant_sites_proportion(partition,
                                                       treeinfo->param_indices[i][j],
                                                       x[part][1]))
            {
              assert(pll_errno);
              return PLL_FAILURE;
            }
          }
          break;
        case PLLMOD_OPT_PARAM_FREE_RATES:
          /* update rate categories */
          memcpy(partition->rates, x[part], partition->rate_cats*sizeof(double));
          break;
        case PLLMOD_OPT_PARAM_RATE_WEIGHTS:
        {
          unsigned int highest_weight_state = fixed_var_index[part];
          unsigned int n_weights            = partition->rate_cats;
          double sum_ratios = 1.0;

          double * weights = partition->rate_weights;
          unsigned int i, cur_weight;

          for (i = 0; i < (n_weights - 1); ++i)
            sum_ratios += x[part][i];

          cur_weight = 0;
          for (i = 0; i < (n_weights); ++i)
            if (i != highest_weight_state)
            {
              weights[i] = x[part][cur_weight++] / sum_ratios;
            }
          weights[highest_weight_state] = 1.0 / sum_ratios;
          break;
        }
        default:
          assert(0);
      }

      part++;
    }
  }

  /* compute negative score */
  if(x)
    score = -1 * pllmod_treeinfo_compute_loglh(treeinfo, 0);

  /* copy per-partition likelihood to the output array */
  if (fx)
  {
    j = 0;
    for (i = 0; i < treeinfo->partition_count; ++i)
    {
      if ((treeinfo->params_to_optimize[i] & params_to_optimize) == params_to_optimize)
        fx[j++] = -1 * treeinfo->partition_loglh[i];
    }
  }

  if (converged)
  {
    /* check if there is at least one unconverged partition in *any* thread */
    if (treeinfo->parallel_reduce_cb)
    {
      treeinfo->parallel_reduce_cb(treeinfo->parallel_context,
                                   &unconverged_flag, 1, PLLMOD_COMMON_REDUCE_SUM);
    }
    converged[num_parts] = unconverged_flag > 0. ? 0 : 1;
  }

  return score;
}

double target_subst_params_func_multi(void * p, double ** x, double * fx,
                                      int * converged)
{
  struct treeinfo_opt_params * params = (struct treeinfo_opt_params *) p;

  pllmod_treeinfo_t * treeinfo      = params->treeinfo;
  unsigned int num_parts            = params->num_opt_partitions;
  unsigned int params_index         = params->params_index;
  unsigned int * subst_free_params  = params->num_free_params;

  double score = -INFINITY;

  /* any partitions which have not converged yet? */
  double unconverged_flag = 0.;

  size_t i, j;
  size_t part = 0;
  for (i = 0; i < treeinfo->partition_count; ++i)
  {
    pll_partition_t * partition = treeinfo->partitions[i];

    if (treeinfo->params_to_optimize[i] & PLLMOD_OPT_PARAM_SUBST_RATES)
    {
      if (!partition || (converged && converged[part]))
      {
        /* partitions has converged, skip it */
        part++;
        continue;
      }

      unconverged_flag = 1.;

      /* function was called solely to check convergence -> no LH computation */
      if (!x)
      {
        part++;
        continue;
      }

      int * symmetries                = treeinfo->subst_matrix_symmetries[i];
      unsigned int states             = partition->states;
      unsigned int subst_params       = (states * (states-1))/2;
      double *subst_rates             = partition->subst_params[params_index];

      /* update subst rates */
      if (symmetries)
      {
        /* assign values to the substitution rates */
        size_t l, k = 0;
        for (l = 0; l <= subst_free_params[part]; ++l)
        {
          double next_value =
                   (l == (unsigned int)symmetries[subst_params - 1]) ? 1.0 : x[part][k++];
          for (j = 0; j < subst_params; j++)
          {
            if ((unsigned int)symmetries[j] == l)
            {
              subst_rates[j] = next_value;
            }
          }
        }
      }
      else
      {
        memcpy (subst_rates, x[part], ((size_t)subst_params - 1) * sizeof(double));
      }

      /* important!! invalidate eigen-decomposition */
      partition->eigen_decomp_valid[params_index] = 0;

      part++;
    }
  }

  /* compute negative score */
  if(x)
    score = -1 * pllmod_treeinfo_compute_loglh(treeinfo, 0);

  /* copy per-partition likelihood to the output array */
  if (fx)
  {
    j = 0;
    for (i = 0; i < treeinfo->partition_count; ++i)
      if (treeinfo->params_to_optimize[i] & PLLMOD_OPT_PARAM_SUBST_RATES)
        fx[j++] = -1 * treeinfo->partition_loglh[i];
  }

  if (converged)
  {
    /* check if there is at least one unconverged partition in *any* thread */
    if (treeinfo->parallel_reduce_cb)
    {
      treeinfo->parallel_reduce_cb(treeinfo->parallel_context,
                                   &unconverged_flag, 1, PLLMOD_COMMON_REDUCE_SUM);
    }
    converged[num_parts] = unconverged_flag > 0. ? 0 : 1;
  }

  return score;
}


double target_freqs_func_multi(void * p, double ** x, double * fx,
                               int * converged)
{
  struct treeinfo_opt_params * params = (struct treeinfo_opt_params *) p;

  pllmod_treeinfo_t * treeinfo      = params->treeinfo;
  unsigned int num_parts            = params->num_opt_partitions;
  unsigned int params_index         = params->params_index;
  unsigned int * highest_freq_state = params->fixed_var_index;

  double score = -INFINITY;

  /* any partitions which have not converged yet? */
  double unconverged_flag = 0.;

  size_t i, j;
  size_t part = 0;
  for (i = 0; i < treeinfo->partition_count; ++i)
  {
    pll_partition_t * partition = treeinfo->partitions[i];

    if (treeinfo->params_to_optimize[i] & PLLMOD_OPT_PARAM_FREQUENCIES)
    {
      if (!partition || (converged && converged[part]))
      {
        /* partitions has converged, skip it */
        part++;
        continue;
      }

      unconverged_flag = 1.;

      /* function was called solely to check convergence -> no LH computation */
      if (!x)
      {
        part++;
        continue;
      }

      unsigned int states             = partition->states;
      double * freqs                  = partition->frequencies[params_index];

      double sum_ratios = 1.0;
      unsigned int cur_index;

      /* update frequencies */
      for (j = 0; j < (states - 1); ++j)
      {
        assert(x[part][j] == x[part][j]);
        sum_ratios += x[part][j];
      }
      cur_index = 0;
      for (j = 0; j < states; ++j)
      {
        if (j != highest_freq_state[part])
        {
          freqs[j] = x[part][cur_index] / sum_ratios;
          cur_index++;
        }
      }
      freqs[highest_freq_state[part]] = 1.0 / sum_ratios;

//      printf("freqs: %f %f %f %f ", freqs[0], freqs[1], freqs[2], freqs[3]);

      /* important!! invalidate eigen-decomposition */
      partition->eigen_decomp_valid[params_index] = 0;

      part++;
    }
  }

  /* compute negative score */
  if (x)
    score = -1 * pllmod_treeinfo_compute_loglh(treeinfo, 0);

  /* copy per-partition likelihood to the output array */
  if (fx)
  {
    j = 0;
    for (i = 0; i < treeinfo->partition_count; ++i)
      if (treeinfo->params_to_optimize[i] & PLLMOD_OPT_PARAM_FREQUENCIES)
        fx[j++] = -1 * treeinfo->partition_loglh[i];
  }

  if (converged)
  {
    /* check if there is at least one unconverged partition in *any* thread */
    if (treeinfo->parallel_reduce_cb)
    {
      treeinfo->parallel_reduce_cb(treeinfo->parallel_context,
                                   &unconverged_flag, 1, PLLMOD_COMMON_REDUCE_SUM);
    }
    converged[num_parts] = unconverged_flag > 0. ? 0 : 1;
  }

  return score;
}
