#include "algo_callback.h"

double target_freqs_func(void *p, double *x)
{
  struct freqs_params * params = (struct freqs_params *) p;

  pll_partition_t * partition     = params->partition;
  pll_utree_t * tree              = params->tree;
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
                                       tree,
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
  pll_utree_t * tree              = params->tree;
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
                                         tree,
                                         params_indices,
                                         1,   /* update pmatrices */
                                         1);  /* update partials */
  return score;
}

double target_alpha_func(void *p, double x)
{
  struct default_params * params = (struct default_params *) p;
  pll_partition_t * partition   = params->partition;
  pll_utree_t * tree            = params->tree;
  unsigned int * params_indices = params->params_indices;

  /* update rate categories */
  if (!pll_compute_gamma_cats (x,
                               partition->rate_cats,
                               partition->rates))
  {
    return PLL_FAILURE;
  }

  /* compute negative score */
  double score = -1 *
                 pllmod_utree_compute_lk(partition,
                                         tree,
                                         params_indices,
                                         1,   /* update pmatrices */
                                         1);  /* update partials */
  return score;
}

double target_pinv_func(void *p, double x)
{
  struct default_params * params = (struct default_params *) p;
  pll_partition_t * partition   = params->partition;
  pll_utree_t * tree            = params->tree;
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
                                         tree,
                                         params_indices,
                                         1,   /* update pmatrices */
                                         1);  /* update partials */
  return score;
}

double target_rates_func(void *p, double *x)
{
  struct rate_weights_params * params = (struct rate_weights_params *) p;
  pll_partition_t * partition   = params->partition;
  pll_utree_t * tree            = params->tree;
  unsigned int * params_indices = params->params_indices;

  /* update rate categories */
  memcpy(partition->rates, x, partition->rate_cats*sizeof(double));

  /* compute negative score */
  double score = -1 *
                 pllmod_utree_compute_lk(partition,
                                         tree,
                                         params_indices,
                                         1,   /* update pmatrices */
                                         1);  /* update partials */
  return score;
}

double target_weights_func(void *p, double *x)
{
  struct rate_weights_params * params = (struct rate_weights_params *) p;

  pll_partition_t * partition       = params->partition;
  pll_utree_t * tree                = params->tree;
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
                                       tree,
                                       params_indices,
                                       0,   /* update pmatrices */
                                       0);  /* update partials */

  return score;
}
