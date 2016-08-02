#include "pllmod_algorithm.h"
#include "algo_callback.h"

static void fill_rates   (double *rates,
                          double *x,
                          int *bt, double *lb, double *ub,
                          double min_rate, double max_rate,
                          unsigned int n_rates);

static void fill_weights (double *weights,
                          unsigned int * highest_weight_index,
                          double *x, int *bt, double *lb, double *ub,
                          unsigned int n_weights);

PLL_EXPORT double pllmod_algo_opt_rates_weights (pll_partition_t * partition,
                                                 pll_utree_t * tree,
                                                 unsigned int * params_indices,
                                                 double min_rate,
                                                 double max_rate,
                                                 double tolerance)
{
  double cur_logl, prev_logl;
  double sum_weightrates, rate_scaler;
  double x[4], lb[4],
        ub[4];
  int bt[4];
  unsigned int i;
  struct rate_weights_params opt_params;
  opt_params.partition = partition;
  opt_params.tree = tree;
  opt_params.params_indices = params_indices;


  /* 2 step BFGS */

  cur_logl = 0;
  do
  {
    prev_logl = cur_logl;

    /* optimize mixture weights */

    double * weights = partition->rate_weights;
    fill_weights(weights, &(opt_params.highest_weight_state), x,
                      bt, lb, ub, partition->rate_cats);

    cur_logl = 1
        * pll_minimize_lbfgsb (x, lb, ub, bt, partition->rate_cats-1,
                               PLLMOD_ALGO_BFGS_FACTR, tolerance,
                               (void *) &opt_params, target_weights_func);


    /* optimize mixture rates */

    fill_rates (partition->rates,
                x, bt, lb, ub,
                min_rate, max_rate,
                partition->rate_cats);

    cur_logl = pll_minimize_lbfgsb(x, lb, ub, bt,
                                   partition->rate_cats,
                                   PLLMOD_ALGO_BFGS_FACTR, tolerance,
                                   (void *) &opt_params, target_rates_func);

  } while (!prev_logl || prev_logl - cur_logl > tolerance);

  /* force constraint sum(weights x rates) = 1.0 */
  sum_weightrates = 0.0;
  for (i=0; i<partition->rate_cats; i++)
    sum_weightrates += partition->rates[i] * partition->rate_weights[i];
  rate_scaler = 1.0 / sum_weightrates;

  for (i=0; i<partition->rate_cats; i++)
      partition->rates[i] *= rate_scaler;

  /* scale branch lengths such that likelihood is conserved */
  pll_utree_scale_branches(tree, sum_weightrates);

  /* update pmatrices and partials according to the new branches */
  cur_logl = -1 *
             pll_utree_compute_lk(partition,
                                  tree,
                                  params_indices,
                                  1,   /* update pmatrices */
                                  1);  /* update partials */

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

  for (i = 0; i < n_rates; i++)
  {
    bt[i] = PLL_LBFGSB_BOUND_BOTH;
    lb[i] = min_rate;
    ub[i] = max_rate;

    if (rates[i] < min_rate || rates[i] > max_rate)
      x[i] = (max_rate + min_rate)/2;
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
  for (i = 1; i < n_weights; i++)
    if (weights[i] > weights[*highest_weight_index])
      *highest_weight_index = i;

  for (i = 0; i < n_weights; i++)
  {
    if (i != *highest_weight_index)
    {
      bt[cur_index] = PLL_LBFGSB_BOUND_BOTH;

      double r = weights[i] / weights[*highest_weight_index];
      lb[cur_index] = PLLMOD_ALGO_MIN_WEIGHT_RATIO;
      ub[cur_index] = PLLMOD_ALGO_MAX_WEIGHT_RATIO;
      x[cur_index] = (r >= lb[i] && r <= ub[i]) ? r : 1.0;
      cur_index++;
    }
  }
}
