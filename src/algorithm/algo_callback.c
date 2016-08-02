#include "algo_callback.h"

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
                 pll_utree_compute_lk(partition,
                                      tree,
                                      params_indices,
                                      1,   /* update pmatrices */
                                      1);  /* update partials */
  return score;
}

double target_weights_func(void *p, double *x)
{
  struct rate_weights_params * params = (struct rate_weights_params *) p;
  pll_partition_t * partition   = params->partition;
  pll_utree_t * tree            = params->tree;
  unsigned int * params_indices = params->params_indices;
  unsigned int highest_weight_state = params->highest_weight_state;
  double sum_ratios = 1.0;
  unsigned int n_weights = partition->rate_cats;
  double weights[PLLMOD_ALGO_MAX_WEIGHTS];
  unsigned int i, cur_weight;

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
  memcpy(partition->rate_weights, weights, partition->rate_cats*sizeof(double));

  /* compute negative score */
  double score = -1 *
                 pll_utree_compute_lk(partition,
                                      tree,
                                      params_indices,
                                      0,   /* update pmatrices */
                                      0);  /* update partials */

  return score;
}
