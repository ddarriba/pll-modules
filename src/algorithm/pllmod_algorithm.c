#include "pllmod_algorithm.h"

#define MAX_WEIGHTS 20

static void fill_rates   (double *rates,
                          double *x,
                          int *bt, double *lb, double *ub,
                          unsigned int n_rates)
{
  unsigned int i;
  for (i = 0; i < n_rates; i++)
  {
    bt[i] = PLL_LBFGSB_BOUND_BOTH;
    lb[i] = 0.01;
    ub[i] = 100;
    assert(rates[i] > lb[i] && rates[i] < ub[i]);
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
      lb[cur_index] = 0.1;
      ub[cur_index] = 10;
      x[cur_index] = (r > lb[i] && r < ub[i]) ? r : 1.0;
      cur_index++;
    }
  }
}

struct rate_weights_params {
  pll_partition_t * partition;
  pll_utree_t * tree;
  unsigned int * params_indices;
  unsigned int highest_weight_state;
};

static double target_rates_func(void *p, double *x)
{
  struct rate_weights_params * params = (struct rate_weights_params *) p;
  pll_partition_t * partition   = params->partition;
  pll_utree_t * tree            = params->tree;
  unsigned int * params_indices = params->params_indices;
  unsigned int i;

  /* update rate categories */
  for (i=0; i<partition->rate_cats; i++)
     partition->rates[i] = x[i];

  /* compute negative score */
  double score = -1 *
                 pll_utree_compute_lk(partition,
                                      tree,
                                      params_indices,
                                      1,   /* update pmatrices */
                                      1);  /* update partials */
  return score;
}

static double target_weights_func(void *p, double *x)
{
  struct rate_weights_params * params = (struct rate_weights_params *) p;
  pll_partition_t * partition   = params->partition;
  pll_utree_t * tree            = params->tree;
  unsigned int * params_indices = params->params_indices;
  unsigned int highest_weight_state = params->highest_weight_state;
  double sum_ratios = 0, sum_weights = 0;
  unsigned int n_weights = partition->rate_cats;
  double weights[MAX_WEIGHTS];
  unsigned int i, cur_weight;

  for (i = 0; i < (n_weights - 1); ++i)
    sum_ratios += x[i];

  cur_weight = 0;
  for (i = 0; i < (n_weights); ++i)
    if (i != highest_weight_state)
    {
      weights[i] = x[cur_weight++] / sum_ratios;
      sum_weights += weights[i];
    }
  weights[highest_weight_state] = 1.0 / sum_ratios;
  sum_weights += weights[highest_weight_state];

  /* update rate categories */
  for (i=0; i<(partition->rate_cats); i++)
     partition->rate_weights[i] = weights[i]/sum_weights;

  /* compute negative score */
  double score = -1 *
                 pll_utree_compute_lk(partition,
                                      tree,
                                      params_indices,
                                      0,   /* update pmatrices */
                                      0);  /* update partials */

  return score;
}

static void scale_branch_length_recursive (pll_utree_t * tree,
                                           double factor)
{
  if (tree)
  {
    tree->length *= factor;
    tree->back->length *= factor;

    if (tree->next)
    {
      scale_branch_length_recursive (tree->next->back, factor);
      scale_branch_length_recursive (tree->next->next->back, factor);
    }
  }
}

PLL_EXPORT double pllmod_algo_opt_rates_weights (pll_partition_t * partition,
                                                 pll_utree_t * tree,
                                                 unsigned int * params_indices,
                                                 double tolerance)
{
  double cur_logl, prev_logl;
  double sumWR, rate_scaler;
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
                               1e9, tolerance,
                               (void *) &opt_params, target_weights_func);


    /* optimize mixture rates */

    fill_rates (partition->rates, x, bt, lb, ub, partition->rate_cats);

    cur_logl = pll_minimize_lbfgsb(x, lb, ub, bt,
                                   partition->rate_cats,
                                   1e9, tolerance,
                                   (void *) &opt_params, target_rates_func);

  } while (!prev_logl || prev_logl - cur_logl > tolerance);

  /* calculate scaler */
  sumWR = 0.0;
  for (i=0; i<partition->rate_cats; i++)
    sumWR += partition->rates[i] * partition->rate_weights[i];
  rate_scaler = 1.0 / sumWR;

  for (i=0; i<partition->rate_cats; i++)
      partition->rates[i] *= rate_scaler;

  /* scale branch lengths */
  scale_branch_length_recursive(tree, sumWR);
  if (tree->back->next)
  {
    scale_branch_length_recursive(tree->back->next->back, sumWR);
    scale_branch_length_recursive(tree->back->next->next->back, sumWR);
  }


   cur_logl = -1 *
              pll_utree_compute_lk(partition,
                                   tree,
                                   params_indices,
                                   1,   /* update pmatrices */
                                   1);  /* update partials */

  return cur_logl;
}
