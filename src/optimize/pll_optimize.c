/*
 Copyright (C) 2015-18 Diego Darriba, Alexey Kozlov

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

 Contact: Alexey Kozlov <Alexey.Kozlov@h-its.org>,
 Exelixis Lab, Heidelberg Instutute for Theoretical Studies
 Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
 */

 /**
  * @file pll_optimize.c
  *
  * @brief Optimization algorithms
  *
  * This file implements high-level optimization algorithms. Most of the
  * functions listed here provide wrappers to call low-level algorithms in
  * `opt_algorithms.c` given a complex structure containing PLL features such
  * as partitions (`pll_partition_t`) or trees (`pll_unode_t`).
  *
  * @author Diego Darriba
  * @author Alexey Kozlov
  */

#include "pll_optimize.h"
#include "lbfgsb/lbfgsb.h"
#include "../pllmod_common.h"

#define BETTER_LL_TRESHOLD 1e-13

/*
 * Note: Compile with flag _ULTRACHECK for checking pre/postconditions
 *       way more thoroughly. This may slow down the execution.
 */
static inline int is_nan(double v)
{
  return v!=v;
}

static int v_int_max (int * v, int n)
{
  int i, max = v[0];
  for (i = 1; i < n; i++)
    if (v[i] > max)
      max = v[i];
  return max;
}

static inline int d_equals(double a, double b)
{
  return (fabs(a-b) < 1e-10);
}

static inline int check_loglh_improvement(int opt_method)
{
  return (opt_method == PLLMOD_OPT_BLO_NEWTON_OLDSAFE ||
          opt_method == PLLMOD_OPT_BLO_NEWTON_SAFE) ? 1 : 0;
}

static int set_x_to_parameters(pll_optimize_options_t * params,
                               double *x)
{
  pll_partition_t * partition = params->lk_params.partition;
  pll_operation_t * operations = params->lk_params.operations;
  double * branch_lengths = params->lk_params.branch_lengths;
  const unsigned int * matrix_indices = params->lk_params.matrix_indices;
  unsigned int params_index = params->params_index;
  const unsigned int * params_indices = params->lk_params.params_indices;
  unsigned int n_branches, n_inner_nodes;
  double * xptr = x;

  if (params->lk_params.rooted)
  {
    n_branches = 2 * partition->tips - 2;
    n_inner_nodes = partition->tips - 1;
  }
  else
  {
    n_branches = 2 * partition->tips - 3;
    n_inner_nodes = partition->tips - 2;
  }

  /* update substitution rate parameters */
  if (params->which_parameters & PLLMOD_OPT_PARAM_SUBST_RATES)
  {
    int * symm;
    int n_subst_rates;
    double * subst_rates;

    symm = params->subst_params_symmetries;
    n_subst_rates = partition->states * (partition->states - 1) / 2;
    if ((subst_rates = (double *) malloc (
        (size_t) n_subst_rates * sizeof(double))) == NULL)
    {
      pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                "Cannot allocate memory for substitution rate parameters");
      return PLL_FAILURE;
    }

    assert(subst_rates);

    if (symm)
    {
      int i, j, k;
      int n_subst_free_params = 0;

      /* compute the number of free parameters */
      n_subst_free_params = v_int_max (symm, n_subst_rates);

      /* assign values to the substitution rates */
      k = 0;
      for (i = 0; i <= n_subst_free_params; i++)
      {
        double next_value = (i == symm[n_subst_rates - 1]) ? 1.0 : xptr[k++];
        for (j = 0; j < n_subst_rates; j++)
          if (symm[j] == i)
          {
            subst_rates[j] = next_value;
          }
      }
      xptr += n_subst_free_params;
    }
    else
    {
      memcpy (subst_rates, xptr, ((size_t)n_subst_rates - 1) * sizeof(double));
      subst_rates[n_subst_rates - 1] = 1.0;
      xptr += n_subst_rates-1;
    }

    pll_set_subst_params (partition,
                          params_index,
                          subst_rates);
    free (subst_rates);
  }

  /* update stationary frequencies */
  if (params->which_parameters & PLLMOD_OPT_PARAM_FREQUENCIES)
  {
    unsigned int i;
    unsigned int n_states = partition->states;
    unsigned int cur_index;
    double sum_ratios = 1.0;
    double *freqs;
    if ((freqs = (double *) malloc ((size_t) n_states * sizeof(double))) == NULL)
    {
      pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                       "Cannot allocate memory for frequencies");
      return PLL_FAILURE;
    }

    for (i = 0; i < (n_states - 1); ++i)
    {
      assert(!is_nan(xptr[i]));
      sum_ratios += xptr[i];
    }
    cur_index = 0;
    for (i = 0; i < (n_states); ++i)
    {
      if (i != params->highest_freq_state)
      {
        freqs[i] = xptr[cur_index] / sum_ratios;
        cur_index++;
      }
    }
    freqs[params->highest_freq_state] = 1.0 / sum_ratios;

    pll_set_frequencies (partition,
                         params_index,
                         freqs);
    free (freqs);
    xptr += (n_states - 1);
  }
  /* update proportion of invariant sites */
  if (params->which_parameters & PLLMOD_OPT_PARAM_PINV)
  {
    assert(!is_nan(xptr[0]));
    unsigned int i;
    for (i = 0; i < (partition->rate_cats); ++i)
    {
      if (!pll_update_invariant_sites_proportion (partition,
                                                  params_indices[i],
                                                  xptr[0]))
      {
        return PLL_FAILURE;
      }
    }
    xptr++;
  }
  /* update gamma shape parameter */
  if (params->which_parameters & PLLMOD_OPT_PARAM_ALPHA)
  {
    assert(!is_nan(xptr[0]));
    /* assign discrete rates */
    double * rate_cats;
    if ((rate_cats = malloc ((size_t) partition->rate_cats * sizeof(double)))
        == NULL)
    {
      pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                "Cannot allocate memory for substitution rate categories");
      return PLL_FAILURE;
    }

    params->lk_params.alpha_value = xptr[0];
    if (!pll_compute_gamma_cats (xptr[0], partition->rate_cats,
                                 rate_cats, PLL_GAMMA_RATES_MEAN))
    {
      return PLL_FAILURE;
    }
    pll_set_category_rates (partition, rate_cats);

    free(rate_cats);
    xptr++;
  }

  /* update free rates */
  if (params->which_parameters & PLLMOD_OPT_PARAM_FREE_RATES)
  {
    pll_set_category_rates (partition, xptr);
    xptr += params->lk_params.partition->rate_cats;
  }

  /* update rate weights */
  if (params->which_parameters & PLLMOD_OPT_PARAM_RATE_WEIGHTS)
  {
    unsigned int i;
    unsigned int rate_cats = params->lk_params.partition->rate_cats;
    unsigned int cur_index;
    double sum_ratios = 1.0;
    double *weights;
    if ((weights = (double *) malloc ((size_t) rate_cats * sizeof(double)))
        == NULL)
    {
      pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                "Cannot allocate memory for substitution rate weights");
      return PLL_FAILURE;
    }

    for (i = 0; i < (rate_cats - 1); ++i)
    {
      assert(!is_nan (xptr[i]));
      sum_ratios += xptr[i];
    }
    cur_index = 0;
    for (i = 0; i < (rate_cats); ++i)
    {
      if (i != params->highest_weight_state)
      {
        weights[i] = xptr[cur_index] / sum_ratios;
        cur_index++;
      }
    }
    weights[params->highest_weight_state] = 1.0 / sum_ratios;
    pll_set_category_weights (partition, weights);
    free (weights);
    xptr += (rate_cats - 1);
  }

  /* update all branch lengths */
  if (params->which_parameters & PLLMOD_OPT_PARAM_BRANCHES_ALL)
  {
    /* assign branch lengths */
    memcpy (branch_lengths, xptr, (size_t)n_branches * sizeof(double));
    xptr += n_branches;
  }

  /* update single branch */
  if (params->which_parameters & PLLMOD_OPT_PARAM_BRANCHES_SINGLE)
   {
    assert(!is_nan(xptr[0]));
     /* assign branch length */
     *branch_lengths = *xptr;
     pll_update_prob_matrices (partition,
                               params_indices,
                               &params->lk_params.where.unrooted_t.edge_pmatrix_index,
                               xptr,
                               1);
     xptr++;
   }
  else
   {
       pll_update_prob_matrices (partition,
                                 params_indices,
                                 matrix_indices,
                                 branch_lengths,
                                 n_branches);

       pll_update_partials (partition, operations, n_inner_nodes);
   }
  return PLL_SUCCESS;
}

static void utree_derivative_func (void * parameters, double proposal,
                                     double *df, double *ddf)
{
  pll_newton_tree_params_t * params = (pll_newton_tree_params_t *) parameters;
  pll_compute_likelihood_derivatives (
      params->partition,
      params->tree->scaler_index,
      params->tree->back->scaler_index,
      proposal,
      params->params_indices,
      params->sumtable, df, ddf);
}


static double compute_negative_lnl_unrooted (void * p, double *x)
{
  pll_optimize_options_t * params = (pll_optimize_options_t *) p;
  pll_partition_t * partition = params->lk_params.partition;
  double score;

  if (x && !set_x_to_parameters(params, x))
    return (double) -INFINITY;

  if (params->lk_params.rooted)
  {
    score = -1
        * pll_compute_root_loglikelihood (
            partition,
            params->lk_params.where.rooted_t.root_clv_index,
            params->lk_params.where.rooted_t.scaler_index,
            params->lk_params.params_indices,
            NULL);
  }
  else
  {
    score = -1
        * pll_compute_edge_loglikelihood (
            partition,
            params->lk_params.where.unrooted_t.parent_clv_index,
            params->lk_params.where.unrooted_t.parent_scaler_index,
            params->lk_params.where.unrooted_t.child_clv_index,
            params->lk_params.where.unrooted_t.child_scaler_index,
            params->lk_params.where.unrooted_t.edge_pmatrix_index,
            params->lk_params.params_indices,
            NULL);
  }

  return score;
} /* compute_lnl_unrooted */

static unsigned int count_n_free_variables (pll_optimize_options_t * params)
{
  unsigned int num_variables = 0;
  pll_partition_t * partition = params->lk_params.partition;

  /* count number of variables for dynamic allocation */
  if (params->which_parameters & PLLMOD_OPT_PARAM_SUBST_RATES)
  {
    int n_subst_rates = partition->states * (partition->states - 1) / 2;
    num_variables +=
        params->subst_params_symmetries ?
            (unsigned int) v_int_max (params->subst_params_symmetries,
                                      n_subst_rates) :
            (unsigned int) n_subst_rates - 1;
  }
  if (params->which_parameters & PLLMOD_OPT_PARAM_FREQUENCIES)
    num_variables += partition->states - 1;
  num_variables += (params->which_parameters & PLLMOD_OPT_PARAM_PINV) != 0;
  num_variables += (params->which_parameters & PLLMOD_OPT_PARAM_ALPHA) != 0;
  if (params->which_parameters & PLLMOD_OPT_PARAM_FREE_RATES)
    num_variables += partition->rate_cats;
  if (params->which_parameters & PLLMOD_OPT_PARAM_RATE_WEIGHTS)
    num_variables += partition->rate_cats - 1;
  num_variables += (params->which_parameters & PLLMOD_OPT_PARAM_BRANCHES_SINGLE)
      != 0;
  if (params->which_parameters & PLLMOD_OPT_PARAM_BRANCHES_ALL)
  {
    unsigned int num_branch_lengths =
        params->lk_params.rooted ?
            (2 * partition->tips - 3) : (2 * partition->tips - 2);
    num_variables += num_branch_lengths;
  }
  return num_variables;
} /* count_n_free_variables */


/******************************************************************************/
/* BRENT'S OPTIMIZATION */
/******************************************************************************/

static double brent_target(void * p, double x)
{
  double score = compute_negative_lnl_unrooted(p, &x);
  return score;
}

/**
 * Optimize one dimension variable with Brent algorithm within a defined range.
 * Target function minimizes the negative likelihood score (i.e., a double
 * precision positive value) given the input parameters.
 * The optimal parameter value is updated in `params`.
 *
 * @param[in,out]  params optimization parameters structure
 * @param      umin   lower bound for target variable
 * @param      umax   upper bound for target variable
 *
 * @return    the negative likelihood score
 */
PLL_EXPORT double pllmod_opt_optimize_onedim(pll_optimize_options_t * params,
                                             double umin,
                                             double umax)
{
  double score = 0;

  /* Brent parameters */
  double xmin;
  double xguess;
  double xmax;
  double f2x;

  switch (params->which_parameters)
  {
    case PLLMOD_OPT_PARAM_ALPHA:
      xguess = params->lk_params.alpha_value;
      xmin   = (umin>0)?umin:PLLMOD_OPT_MIN_ALPHA;
      xmax   = (umax>0)?umax:PLLMOD_OPT_MAX_ALPHA;
      break;
    case PLLMOD_OPT_PARAM_PINV:
      xguess = params->lk_params.partition->prop_invar[params->params_index];
      xmin   = (umin>0)?umin:PLLMOD_OPT_MIN_PINV;
      xmax   = (umax>0)?umax:PLLMOD_OPT_MAX_PINV;
      break;
    case PLLMOD_OPT_PARAM_BRANCHES_SINGLE:
      xguess = params->lk_params.branch_lengths[0];
      xmin   = (umin>0)?umin:PLLMOD_OPT_MIN_BRANCH_LEN;
      xmax   = (umax>0)?umax:PLLMOD_OPT_MAX_BRANCH_LEN;
      break;
    default:
      /* unavailable or multiple parameter */
      return (double) -INFINITY;
  }

  double xres = pllmod_opt_minimize_brent(xmin, xguess, xmax,
                                          params->pgtol,
                                          &score,
                                          &f2x,
                                          (void *) params,
                                          &brent_target);
  set_x_to_parameters(params, &xres);

  return score;
} /* pll_optimize_parameters_onedim */

/******************************************************************************/
/* L-BFGS-B OPTIMIZATION */
/******************************************************************************/

/**
 * Optimize multi-dimensional variable with L-BFGS-B algorithm within a defined
 * range.
 * Target function minimizes the negative likelihood score (i.e., a double
 * precision positive value) given the input parameters.
 * The optimal parameter values are updated in `params`.
 *
 * @param[in,out]  params optimization parameters structure
 * @param  umin   array containing lower bounds for target variables
 * @param  umax   array containing upper bounds for target variables
 *
 * @return        the negative likelihood score
 */
PLL_EXPORT double pllmod_opt_optimize_multidim (pll_optimize_options_t * params,
                                                double *umin,
                                                double *umax)
{
  unsigned int i;
  pll_partition_t * partition = params->lk_params.partition;

  /* L-BFGS-B parameters */
  //  double initial_score;
  unsigned int num_variables;
  double score = 0;
  double *x, *lower_bounds, *upper_bounds;
  int *bound_type;

  /* ensure that the 2 branch optimization modes are not set together */
  assert(!((params->which_parameters & PLLMOD_OPT_PARAM_BRANCHES_ALL)
      && (params->which_parameters & PLLMOD_OPT_PARAM_BRANCHES_SINGLE)));

  num_variables = count_n_free_variables (params);

  x = (double *) calloc ((size_t) num_variables, sizeof(double));
  lower_bounds = (double *) calloc ((size_t) num_variables, sizeof(double));
  upper_bounds = (double *) calloc ((size_t) num_variables, sizeof(double));
  bound_type = (int *) calloc ((size_t) num_variables, sizeof(int));

  if (!(x && lower_bounds && upper_bounds && bound_type))
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
              "Cannot allocate memory for l-bfgs-b parameters");
    if (x)
      free (x);
    if (lower_bounds)
      free (lower_bounds);
    if (upper_bounds)
      free (upper_bounds);
    if (bound_type)
      free (bound_type);
    return (double) -INFINITY;
  }

  {
    int * nbd_ptr = bound_type;
    /* effective boundaries */
    double * l_ptr = lower_bounds, *u_ptr = upper_bounds;
    /* user defined boundaries */
    double * ul_ptr = umin, *uu_ptr = umax;
    unsigned int check_n = 0;

    /* substitution rate parameters */
    if (params->which_parameters & PLLMOD_OPT_PARAM_SUBST_RATES)
    {
      unsigned int n_subst_rates;
      unsigned int n_subst_free_params;

      n_subst_rates = partition->states * (partition->states - 1) / 2;
      if (params->subst_params_symmetries)
      {
        n_subst_free_params =(unsigned int) v_int_max (
                                         params->subst_params_symmetries,
                                         (int) n_subst_rates);
      }
      else
      {
        n_subst_free_params = n_subst_rates - 1;
      }

      int current_rate = 0;
      for (i = 0; i < n_subst_free_params; i++)
      {
        nbd_ptr[i] = PLLMOD_OPT_LBFGSB_BOUND_BOTH;
        unsigned int j = i;
        if (params->subst_params_symmetries)
        {
          if (params->subst_params_symmetries[n_subst_rates-1] == current_rate)
            current_rate++;
          for (j=0; j<n_subst_rates; j++)
          {
            if (params->subst_params_symmetries[j] == current_rate)
              break;
          }
          current_rate++;
        }

        x[check_n + i] = partition->subst_params[params->params_index][j];
        l_ptr[i] = ul_ptr?(*(ul_ptr++)):PLLMOD_OPT_MIN_SUBST_RATE;
        u_ptr[i] = uu_ptr?(*(uu_ptr++)):PLLMOD_OPT_MAX_SUBST_RATE;
      }
      nbd_ptr += n_subst_free_params;
      l_ptr += n_subst_free_params;
      u_ptr += n_subst_free_params;
      check_n += n_subst_free_params;
    }

    /* stationary frequency parameters */
    if (params->which_parameters & PLLMOD_OPT_PARAM_FREQUENCIES)
    {
      unsigned int states = params->lk_params.partition->states;
      unsigned int n_freqs_free_params = states - 1;
      unsigned int cur_index;

      double * frequencies =
          params->lk_params.partition->frequencies[params->params_index];

      params->highest_freq_state = 3;
      for (i = 1; i < states; i++)
              if (frequencies[i] > frequencies[params->highest_freq_state])
                params->highest_freq_state = i;

      cur_index = 0;
      for (i = 0; i < states; i++)
      {
        if (i != params->highest_freq_state)
        {
          nbd_ptr[cur_index] = PLLMOD_OPT_LBFGSB_BOUND_BOTH;
          x[check_n + cur_index] = frequencies[i]
              / frequencies[params->highest_freq_state];
          l_ptr[cur_index] = ul_ptr?(*(ul_ptr++)):PLLMOD_OPT_MIN_FREQ;
          u_ptr[cur_index] = uu_ptr?(*(uu_ptr++)):PLLMOD_OPT_MAX_FREQ;
          cur_index++;
        }
      }
      check_n += n_freqs_free_params;
      nbd_ptr += n_freqs_free_params;
      l_ptr += n_freqs_free_params;
      u_ptr += n_freqs_free_params;
    }

    /* proportion of invariant sites */
    if (params->which_parameters & PLLMOD_OPT_PARAM_PINV)
    {
      *nbd_ptr = PLLMOD_OPT_LBFGSB_BOUND_BOTH;
      x[check_n] = partition->prop_invar[params->params_index];
      *l_ptr = ul_ptr?(*(ul_ptr++)):PLLMOD_OPT_MIN_PINV + PLL_LBFGSB_ERROR;
      *u_ptr = uu_ptr?(*(uu_ptr++)):PLLMOD_OPT_MAX_PINV;
      check_n++;
      nbd_ptr++;
      l_ptr++;
      u_ptr++;
    }

    /* gamma shape parameter */
    if (params->which_parameters & PLLMOD_OPT_PARAM_ALPHA)
    {
      *nbd_ptr = PLLMOD_OPT_LBFGSB_BOUND_BOTH;
      x[check_n] = params->lk_params.alpha_value;
      *l_ptr = ul_ptr?(*(ul_ptr++)):PLLMOD_OPT_MIN_ALPHA;
      *u_ptr = uu_ptr?(*(uu_ptr++)):PLLMOD_OPT_MAX_ALPHA;
      check_n++;
      nbd_ptr++;
      l_ptr++;
      u_ptr++;
    }

    /* update free rates */
      if (params->which_parameters & PLLMOD_OPT_PARAM_FREE_RATES)
      {
        unsigned int n_cats = params->lk_params.partition->rate_cats;
        for (i=0; i<n_cats; i++)
        {
          x[check_n + i]  = params->lk_params.partition->rates[i];
          l_ptr[i] = ul_ptr?(*(ul_ptr++)):PLLMOD_OPT_MIN_RATE;
          u_ptr[i] = uu_ptr?(*(uu_ptr++)):PLLMOD_OPT_MAX_RATE;
          nbd_ptr[i] = PLLMOD_OPT_LBFGSB_BOUND_BOTH;
        }
        check_n += n_cats;
        nbd_ptr += (int) n_cats;
        l_ptr   += (int) n_cats;
        u_ptr   += (int) n_cats;
      }

      if (params->which_parameters & PLLMOD_OPT_PARAM_RATE_WEIGHTS)
      {
      unsigned int rate_cats = params->lk_params.partition->rate_cats;
      unsigned int n_weights_free_params = rate_cats - 1;
      unsigned int cur_index;

      double * rate_weights = params->lk_params.partition->rate_weights;

      params->highest_weight_state = rate_cats - 1;
      for (i = 1; i < rate_cats; i++)
        if (rate_weights[i] > rate_weights[params->highest_weight_state])
          params->highest_weight_state = i;

      cur_index = 0;
      for (i = 0; i < rate_cats; i++)
      {
        if (i != params->highest_weight_state)
        {
          nbd_ptr[cur_index] = PLLMOD_OPT_LBFGSB_BOUND_BOTH;
          x[check_n + cur_index] = rate_weights[i]
              / rate_weights[params->highest_weight_state];
          l_ptr[cur_index] = ul_ptr ? (*(ul_ptr++)) : PLLMOD_OPT_MIN_RATE_WEIGHT;
          u_ptr[cur_index] = uu_ptr ? (*(uu_ptr++)) : PLLMOD_OPT_MAX_RATE_WEIGHT;
          cur_index++;
        }
      }
      check_n += n_weights_free_params;
      nbd_ptr += n_weights_free_params;
      l_ptr += n_weights_free_params;
      u_ptr += n_weights_free_params;
    }

    /* topology (UNIMPLEMENTED) */
    if (params->which_parameters & PLLMOD_OPT_PARAM_TOPOLOGY)
    {
      return PLL_FAILURE;
    }

    /* single branch length */
    if (params->which_parameters & PLLMOD_OPT_PARAM_BRANCHES_SINGLE)
    {
        nbd_ptr [check_n]= PLLMOD_OPT_LBFGSB_BOUND_LOWER;
        x[check_n] = params->lk_params.branch_lengths[0];
        l_ptr[check_n] = ul_ptr?(*(ul_ptr++)):PLLMOD_OPT_MIN_BRANCH_LEN;
        u_ptr[check_n] = uu_ptr?(*(uu_ptr++)):PLLMOD_OPT_MAX_BRANCH_LEN;
        check_n++;
        nbd_ptr++;
        l_ptr++;
        u_ptr++;
    }

    /* all branches */
    if (params->which_parameters & PLLMOD_OPT_PARAM_BRANCHES_ALL)
    {
      unsigned int num_branch_lengths =
          params->lk_params.rooted ?
              (2 * partition->tips - 3) : (2 * partition->tips - 2);
      for (i = 0; i < num_branch_lengths; i++)
      {
        nbd_ptr[i] = PLLMOD_OPT_LBFGSB_BOUND_LOWER;
        x[check_n + i] = params->lk_params.branch_lengths[i];
        l_ptr[check_n + i] = ul_ptr?(*(ul_ptr++)):PLLMOD_OPT_MIN_BRANCH_LEN;
        u_ptr[check_n + i] = uu_ptr?(*(uu_ptr++)):PLLMOD_OPT_MAX_BRANCH_LEN;
      }
      check_n += num_branch_lengths;
      nbd_ptr += num_branch_lengths;
      l_ptr += num_branch_lengths;
      u_ptr += num_branch_lengths;
    }
    assert(check_n == num_variables);
  }

  score = pllmod_opt_minimize_lbfgsb(x, lower_bounds, upper_bounds, bound_type,
                              num_variables, params->factr, params->pgtol,
                              params, compute_negative_lnl_unrooted);

  free (x);
  free (lower_bounds);
  free (upper_bounds);
  free (bound_type);

  if (is_nan(score))
  {
    score = (double) -INFINITY;
    if (!pll_errno)
    {
      pllmod_set_error(PLLMOD_OPT_ERROR_LBFGSB_UNKNOWN,
                       "Unknown LBFGSB error");
    }
  }

  return score;
} /* pllmod_opt_optimize_multidim */

/******************************************************************************/
/* GENERIC */
/******************************************************************************/

static void update_partials_and_scalers(pll_partition_t ** partitions,
                                        size_t partition_count,
                                        pll_unode_t * parent,
                                        pll_unode_t * right_child,
                                        pll_unode_t * left_child)
{
  pll_operation_t op;
  size_t p;

  /* set CLV */
  op.parent_clv_index    = parent->clv_index;
  op.parent_scaler_index = parent->scaler_index;
  op.child1_clv_index    = right_child->back->clv_index;
  op.child1_matrix_index = right_child->back->pmatrix_index;
  op.child1_scaler_index = right_child->back->scaler_index;
  op.child2_clv_index    = left_child->back->clv_index;
  op.child2_matrix_index = left_child->back->pmatrix_index;
  op.child2_scaler_index = left_child->back->scaler_index;

  for (p = 0; p < partition_count; ++p)
  {
    /* skip remote partitions */
    if (!partitions[p])
      continue;

    pll_update_partials (partitions[p], &op, 1);
  }
}

/* if keep_update, P-matrices are updated after each branch length opt */
static int recomp_iterative (pll_newton_tree_params_t * params,
                              int radius,
                              double * loglikelihood_score,
                              int keep_update)
{
  pll_unode_t *tr_p, *tr_q, *tr_z;
  double xmin,    /* min branch length */
         xguess,  /* initial guess */
         xmax,    /* max branch length */
         xtol,    /* tolerance */
         xres,    /* optimal found branch length */
         xorig;   /* original branch length before optimization */

  tr_p = params->tree;
  tr_q = params->tree->next;
  tr_z = tr_q ? tr_q->next : NULL;
  xorig = tr_p->length;

  /* check branch length integrity */
  assert(d_equals(tr_p->length, tr_p->back->length));

  /* prepare sumtable for current branch */
  pll_update_sumtable (params->partition,
                       tr_p->clv_index,
                       tr_p->back->clv_index,
                       tr_p->scaler_index,
                       tr_p->back->scaler_index,
                       params->params_indices,
                       params->sumtable);

  /* set N-R parameters */
  xmin = params->branch_length_min;
  xmax = params->branch_length_max;
  xtol = params->tolerance;
  xguess = tr_p->length;
  if (xguess < xmin || xguess > xmax)
    xguess = PLLMOD_OPT_DEFAULT_BRANCH_LEN;

  xres = pllmod_opt_minimize_newton_old(xmin, xguess, xmax, xtol,
                              params->max_newton_iters, params,
                              utree_derivative_func);

  if (pll_errno)
    return PLL_FAILURE;

  /* update branch length in the tree structure */
  tr_p->length = tr_p->back->length = xres;

  if (keep_update && fabs(tr_p->length - xorig) > 1e-10)
  {
    /* update pmatrix for the new branch length */
    pll_update_prob_matrices(params->partition,
                             params->params_indices,
                             &(tr_p->pmatrix_index),
                             &xres,1);

    if (check_loglh_improvement(params->opt_method))
    {
      /* check and compare likelihood */
      double eval_loglikelihood = pll_compute_edge_loglikelihood (params->partition,
                                                                  tr_p->clv_index,
                                                                  tr_p->scaler_index,
                                                                  tr_p->back->clv_index,
                                                                  tr_p->back->scaler_index,
                                                                  tr_p->pmatrix_index,
                                                                  params->params_indices,
                                                                  NULL);

      /* check if the optimal found value improves the likelihood score */
      if (eval_loglikelihood >= *loglikelihood_score)
      {
        /* fix new score */
        *loglikelihood_score = eval_loglikelihood;

        /* update branch length in the tree structure */
        tr_p->length = xres;
        tr_p->back->length = tr_p->length;
      }
      else
      {
        /* reset branch length to original value */
        tr_p->length = tr_p->back->length = xorig;

        pll_update_prob_matrices(params->partition,
                                 params->params_indices,
                                 &(tr_p->pmatrix_index),
                                 &tr_p->length,1);
      }
    }
  }

  DBG(" Optimized branch %3d - %3d (%.6f)\n",
      tr_p->clv_index, tr_p->back->clv_index, tr_p->length);

  /* update children */
  if (radius && tr_q && tr_z)
  {
    /* update children 'Q'
     * CLV at P is recomputed with children P->back and Z->back
     * Scaler is updated by subtracting Q->back and adding P->back
     */
    update_partials_and_scalers(&params->partition,
                                1,
                                tr_q,
                                tr_p,
                                tr_z);

    /* eval */
    pll_newton_tree_params_t params_cpy;
    memcpy(&params_cpy, params, sizeof(pll_newton_tree_params_t));
    params_cpy.tree = tr_q->back;
    if (!recomp_iterative (&params_cpy,
                           radius-1,
                           loglikelihood_score,
                           keep_update))
      return PLL_FAILURE;

    /* update children 'Z'
     * CLV at P is recomputed with children P->back and Q->back
     * Scaler is updated by subtracting Z->back and adding Q->back
     */
    update_partials_and_scalers(&params->partition,
                                1,
                                tr_z,
                                tr_q,
                                tr_p);

   /* eval */
    params_cpy.tree = tr_z->back;
    if (!recomp_iterative (&params_cpy,
                           radius-1,
                           loglikelihood_score,
                           keep_update))
      return PLL_FAILURE;

    /* reset to initial state
     * CLV at P is recomputed with children Q->back and Z->back
     * Scaler is updated by subtracting P->back and adding Z->back
     */
    update_partials_and_scalers(&params->partition,
                                1,
                                tr_p,
                                tr_z,
                                tr_q);
  }

  return PLL_SUCCESS;

} /* recomp_iterative */

/**
 * Optimize branch lengths within a certain radius around the virtual rooted
 * using Newton-Raphson minimization algorithm.
 *
 * There are 2 preconditions for this function:
 *
 * 1. When this function is called, CLVs must be up-to-date towards the virtual
 * root defined by `tree`. This condition cannot be checked, so make sure they
 * are correct.
 *
 * 2. Pmatrix indices must be unique for the involved branches, otherwise there
 * will be side effects when the branches are optimized.
 *
 *
 * `keep_update` determines whether after optimizing a branch, the resulting
 * length is updated in the tree and partition structures before proceeding to
 * the next branch. Otherwise, all branches are optimized given the original
 * branch lengths and CLVs and updated all together at the end of each iteration.
 * In general, `keep_update` provides better fitness, but the results may not be
 * reproducible if several branches are optimized in parallel.
 *
 * @param[in,out]  partition         the PLL partition structure
 * @param[in,out]  tree              the PLL unrotted tree structure
 * @param  params_indices    the indices of the parameter sets
 * @param  branch_length_min lower bound for branch lengths
 * @param  branch_length_max upper bound for branch lengths
 * @param  tolerance         tolerance for Newton-Raphson algorithm
 * @param  smoothings        number of iterations over the branches
 * @param  radius            radius from the virtual root
 * @param  keep_update       if true, branch lengths are iteratively updated in the tree structure
 *
 * @return                   the likelihood score after optimizing branch lengths
 */
PLL_EXPORT double pllmod_opt_optimize_branch_lengths_local (
                                              pll_partition_t * partition,
                                              pll_unode_t * tree,
                                              const unsigned int * params_indices,
                                              double branch_length_min,
                                              double branch_length_max,
                                              double tolerance,
                                              int smoothings,
                                              int radius,
                                              int keep_update)
{
  unsigned int iters;
  double loglikelihood = 0.0, new_loglikelihood;
  unsigned int sites_alloc;

  /*
   * preconditions:
   *    (1) CLVs must be updated towards 'tree'
   *    (2) Pmatrix indices must be **unique** for each branch
   */

  pllmod_reset_error();

  if (radius < PLLMOD_OPT_BRLEN_OPTIMIZE_ALL)
  {
    pllmod_set_error(PLLMOD_OPT_ERROR_NEWTON_BAD_RADIUS,
                     "Invalid radius for branch length optimization");
    return (double)PLL_FAILURE;
  }

  /* get the initial likelihood score */
  loglikelihood = pll_compute_edge_loglikelihood (partition,
                                                  tree->back->clv_index,
                                                  tree->back->scaler_index,
                                                  tree->clv_index,
                                                  tree->scaler_index,
                                                  tree->pmatrix_index,
                                                  params_indices,
                                                  NULL);

  /* set parameters for N-R optimization */
  pll_newton_tree_params_t params;
  params.partition         = partition;
  params.tree              = tree;
  params.params_indices    = params_indices;
  params.branch_length_min = (branch_length_min>0)?
                              branch_length_min:
                              PLLMOD_OPT_MIN_BRANCH_LEN;
  params.branch_length_max = (branch_length_max>0)?
                              branch_length_max:
                              PLLMOD_OPT_MAX_BRANCH_LEN;
  params.tolerance         = (branch_length_min>0)?
                              branch_length_min/10.0:
                              PLLMOD_OPT_TOL_BRANCH_LEN;
  params.sumtable          = 0;
  params.opt_method        = PLLMOD_OPT_BLO_NEWTON_OLDFAST;
  params.max_newton_iters  = 30;

  /* allocate the sumtable */
  sites_alloc = partition->sites;
  if (partition->attributes & PLL_ATTRIB_AB_FLAG)
    sites_alloc += partition->states;

  if ((params.sumtable = (double *) pll_aligned_alloc(
       sites_alloc * partition->rate_cats * partition->states_padded *
       sizeof(double), partition->alignment)) == NULL)
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for bl opt variables");
    return PLL_FAILURE;
  }

  iters = (unsigned int) smoothings;
  while (iters)
  {
    new_loglikelihood = loglikelihood;

    /* iterate on first edge */
    params.tree = tree;
    if (!recomp_iterative (&params, radius, &new_loglikelihood, keep_update))
    {
      loglikelihood = PLL_FAILURE;
      break;
    }

    if (radius)
    {
      /* iterate on second edge */
      params.tree = tree->back;
      if (!recomp_iterative (&params, radius-1, &new_loglikelihood, keep_update))
      {
        loglikelihood = PLL_FAILURE;
        break;
      }
    }
    /* compute likelihood after optimization */
    new_loglikelihood = pll_compute_edge_loglikelihood (partition,
                                                        tree->back->clv_index,
                                                        tree->back->scaler_index,
                                                        tree->clv_index,
                                                        tree->scaler_index,
                                                        tree->pmatrix_index,
                                                        params_indices,
                                                        NULL);

    DBG("pllmod_opt_optimize_branch_lengths_local: iters %d, old: %f, new: %f\n",
        iters, loglikelihood, new_loglikelihood);

    if (new_loglikelihood - loglikelihood > new_loglikelihood * BETTER_LL_TRESHOLD)
    {
      iters --;

      /* check convergence */
      if (fabs (new_loglikelihood - loglikelihood) < tolerance) iters = 0;

      loglikelihood = new_loglikelihood;
    }
    else
    {
      if (check_loglh_improvement(params.opt_method))
        assert(new_loglikelihood - loglikelihood > new_loglikelihood * 1e-14);
      else
      {
        pllmod_set_error(PLLMOD_OPT_ERROR_NEWTON_WORSE_LK,
                         "Local BL opt converged to a worse likelihood score by %f units",
                         new_loglikelihood - loglikelihood);
        loglikelihood = new_loglikelihood;
        break;
      }
    }
  }

  /* deallocate sumtable */
  pll_aligned_free(params.sumtable);

  return -1*loglikelihood;
} /* pllmod_opt_optimize_branch_lengths_local */

/**
 * Optimize branch lengths using Newton-Raphson minimization algorithm.
 *
 * Check `pllmod_opt_optimize_branch_lengths_local` documentation.
 *
 * @param[in,out]  partition         the PLL partition structure
 * @param[in,out]  tree              the PLL unrotted tree structure
 * @param  params_indices    the indices of the parameter sets
 * @param  branch_length_min lower bound for branch lengths
 * @param  branch_length_max upper bound for branch lengths
 * @param  tolerance         tolerance for Newton-Raphson algorithm
 * @param  smoothings        number of iterations over the branches
 * @param  keep_update       if true, branch lengths are iteratively updated in the tree structure
 *
 * @return                   the likelihood score after optimizing branch lengths
 */
PLL_EXPORT double pllmod_opt_optimize_branch_lengths_iterative (
                                              pll_partition_t * partition,
                                              pll_unode_t * tree,
                                              const unsigned int * params_indices,
                                              double branch_length_min,
                                              double branch_length_max,
                                              double tolerance,
                                              int smoothings,
                                              int keep_update)
{
  double loglikelihood;
  loglikelihood = pllmod_opt_optimize_branch_lengths_local (partition,
                                                 tree,
                                                 params_indices,
                                                 branch_length_min,
                                                 branch_length_max,
                                                 tolerance,
                                                 smoothings,
                                                 PLLMOD_OPT_BRLEN_OPTIMIZE_ALL,
                                                 keep_update);
  return loglikelihood;
} /* pllmod_opt_optimize_branch_lengths_iterative */

/**
 * Compute the likelihood function derivatives for a specific branch length
 *
 * @param parameters  `pll_optimize_options_t` structure
 * @param proposal    the branch length where the derivatives are computed
 * @param df[out]         first derivative of the likelihood function
 * @param ddf[out]        second derivative of the likelihood function
 */
PLL_EXPORT void pllmod_opt_derivative_func(void * parameters,
                                          double proposal,
                                          double *df, double *ddf)
{
  pll_optimize_options_t * params = (pll_optimize_options_t *) parameters;
  pll_compute_likelihood_derivatives (
      params->lk_params.partition,
      params->lk_params.where.unrooted_t.parent_scaler_index,
      params->lk_params.where.unrooted_t.child_scaler_index,
      proposal,
      params->lk_params.params_indices,
      params->sumtable,
      df, ddf);
}

/**
 *  multi-partition optimization routines
 */

/**
 * Compute the likelihood score at a given edge on a multiple partition
 *
 * @param  partitions          list of partitions
 * @param  partition_count     number of partitions in `partitions`
 * @param  parent_clv_index    parent clv index
 * @param  parent_scaler_index parent scaler index
 * @param  child_clv_index     child clv index
 * @param  child_scaler_index  child scaler index
 * @param  matrix_index        matrix index of the edge
 * @param  params_indices      the indices of the parameter sets
 * @param  persite_lnl         per-site likelihoods (if 0, they are omitted)
 * @param  parallel_context    context for parallel computation
 * @param  parallel_reduce_cb  callback function for parallel reduction
 *
 * @return                     the likelihood score at the given edge
 */
PLL_EXPORT double pllmod_opt_compute_edge_loglikelihood_multi(
                                              pll_partition_t ** partitions,
                                              size_t partition_count,
                                              unsigned int parent_clv_index,
                                              int parent_scaler_index,
                                              unsigned int child_clv_index,
                                              int child_scaler_index,
                                              unsigned int matrix_index,
                                              unsigned int ** const params_indices,
                                              double * persite_lnl,
                                              void * parallel_context,
                                              void (*parallel_reduce_cb)(void *,
                                                                         double *,
                                                                         size_t,
                                                                         int))
{
  double total_loglh = 0.;

  size_t p;
  for (p = 0; p < partition_count; ++p)
  {
    /* skip remote partitions */
    if (!partitions[p])
      continue;

    total_loglh += pll_compute_edge_loglikelihood(partitions[p],
                                                  parent_clv_index,
                                                  parent_scaler_index,
                                                  child_clv_index,
                                                  child_scaler_index,
                                                  matrix_index,
                                                  params_indices[p],
                                                  NULL);
  }

  if (parallel_reduce_cb)
    parallel_reduce_cb(parallel_context, &total_loglh, 1, PLLMOD_COMMON_REDUCE_SUM);

  return total_loglh;
}

static void utree_derivative_func_multi (void * parameters, double proposal,
                                         double *df, double *ddf)
{
  pll_newton_tree_params_multi_t * params =
                                (pll_newton_tree_params_multi_t *) parameters;
  size_t p;

  *df = *ddf = 0;

  /* simply iterate over partitions and add up the derivatives */
  for (p = 0; p < params->partition_count; ++p)
  {
    /* skip remote partitions */
    if (!params->partitions[p])
      continue;

    double p_df, p_ddf;
    double s = params->brlen_scalers ? params->brlen_scalers[p] : 1.;
    double p_brlen = s * proposal;
    pll_compute_likelihood_derivatives (params->partitions[p],
                                        params->tree->scaler_index,
                                        params->tree->back->scaler_index,
                                        p_brlen,
                                        params->params_indices[p],
                                        params->precomp_buffers[p],
                                        &p_df, &p_ddf);

    /* chain rule! */
    *df += s * p_df;
    *ddf += s * s * p_ddf;
  }

  if (params->parallel_reduce_cb)
  {
    double d[2] = {*df, *ddf};
    params->parallel_reduce_cb(params->parallel_context, d, 2, PLLMOD_COMMON_REDUCE_SUM);
    *df = d[0];
    *ddf = d[1];
  }
}

static void update_prob_matrices(pll_partition_t ** partitions,
                                 size_t partition_count,
                                 unsigned int ** params_indices,
                                 double * brlen_scalers,
                                 pll_unode_t * node)
{
  unsigned int p;
  for (p = 0; p < partition_count; ++p)
  {
    /* skip remote partitions */
    if (!partitions[p])
      continue;

    const double p_brlen = brlen_scalers ?
        brlen_scalers[p] * node->length : node->length;

    pll_update_prob_matrices(partitions[p],
                             params_indices[p],
                             &(node->pmatrix_index),
                             &p_brlen, 1);
  }
}


static int allocate_buffers(pll_newton_tree_params_multi_t * params)
{
  if (!params->precomp_buffers)
  {
    params->precomp_buffers =  (double **) calloc(params->partition_count,
                                                  sizeof(double *));
    if (!params->precomp_buffers)
      return PLL_FAILURE;

    for (unsigned int p = 0; p < params->partition_count; ++p)
    {
      const pll_partition_t * partition = params->partitions[p];

      /* skip remote partitions */
      if (!partition)
        continue;

      unsigned int sites_alloc = partition->sites;
      if (partition->attributes & PLL_ATTRIB_AB_FLAG)
        sites_alloc += partition->states;

      params->precomp_buffers[p] = (double *) pll_aligned_alloc(
                 sites_alloc * partition->rate_cats * partition->states_padded *
                 sizeof(double), partition->alignment);

      if (!params->precomp_buffers[p])
        return PLL_FAILURE;
    }
  }

  if (!params->brlen_buffers && params->opt_method == PLLMOD_OPT_BLO_NEWTON_FALLBACK)
  {
    params->brlen_buffers =  (double **) calloc(params->partition_count,
                                                sizeof(double *));
    if (!params->brlen_buffers)
      return PLL_FAILURE;

    // not very elegant...
    unsigned int branch_count = 0;
    for (size_t p = 0; p < params->partition_count; ++p)
    {
      if (params->partitions[p])
      {
        branch_count = 2 * params->partitions[p]->tips - 3;
        break;
      }
    }
    assert(branch_count);

    if (params->brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED)
    {
      for (size_t p = 0; p < params->partition_count; ++p)
      {
        const pll_partition_t * partition = params->partitions[p];

        /* skip remote partitions */
        if (!partition)
          continue;

        params->brlen_buffers[p] = (double *) calloc(branch_count,
            sizeof(double));

        if (!params->brlen_buffers[p])
          return PLL_FAILURE;
      }
    }
    else
    {

      params->brlen_buffers[0] = (double *) calloc(branch_count,
          sizeof(double));

      if (!params->brlen_buffers[0])
        return PLL_FAILURE;
    }
  }

  return PLL_SUCCESS;
}

/* if keep_update, P-matrices are updated after each branch length opt */
static int recomp_iterative_multi (pll_newton_tree_params_multi_t * params,
                                   int radius,
                                   double * loglikelihood_score,
                                   int keep_update)
{
  pll_unode_t *tr_p, *tr_q, *tr_z;
  size_t p;
  double xmin,    /* min branch length */
         xguess,  /* initial guess */
         xmax,    /* max branch length */
         xtol,    /* tolerance */
         xres,    /* optimal found branch length */
         xorig;   /* original branch length before optimization */

  tr_p = params->tree;
  tr_q = params->tree->next;
  tr_z = tr_q ? tr_q->next : NULL;
  xorig = tr_p->length;

  /* check branch length integrity */
  assert(d_equals(tr_p->length, tr_p->back->length));

  /* prepare sumtable for current branch */
  for (p = 0; p < params->partition_count; ++p)
  {
    /* skip remote partitions */
    if (!params->partitions[p])
      continue;

    pll_update_sumtable (params->partitions[p],
                         tr_p->clv_index,
                         tr_p->back->clv_index,
                         tr_p->scaler_index,
                         tr_p->back->scaler_index,
                         params->params_indices[p],
                         params->precomp_buffers[p]);
  }

  /* set N-R parameters */
  xmin = params->branch_length_min;
  xmax = params->branch_length_max;
  xtol = params->tolerance;
  xguess = tr_p->length;

  switch (params->opt_method)
  {
    case PLLMOD_OPT_BLO_NEWTON_FAST:
    case PLLMOD_OPT_BLO_NEWTON_SAFE:
    {
      xres = pllmod_opt_minimize_newton(xmin, xguess, xmax, xtol,
                                         params->max_newton_iters, params,
                                         utree_derivative_func_multi);
    }
    break;
    case PLLMOD_OPT_BLO_NEWTON_FALLBACK:
      // TODO: adapt for unlinked branches
      params->brlen_buffers[0][tr_p->pmatrix_index] = tr_p->length;
      assert(0);
    case PLLMOD_OPT_BLO_NEWTON_OLDFAST:
    case PLLMOD_OPT_BLO_NEWTON_OLDSAFE:
    {
      xres = pllmod_opt_minimize_newton_old(xmin, xguess, xmax, xtol,
                                         params->max_newton_iters, params,
                                         utree_derivative_func_multi);
    }
    break;
    case PLLMOD_OPT_BLO_NEWTON_GLOBAL:
    {
      assert(0);
    }
    break;
    default:
      assert(0);
  }

  if (pll_errno == PLLMOD_OPT_ERROR_NEWTON_LIMIT)
  {
    /* NR optimization failed to converge:
     * - if LH improvement check is enabled, it is safe to keep
     *   the branch length from the last iteration
     * - otherwise, we must reset branch length to the original value
     *   to avoid getting worse LH score in the end
     * */

    DBG("NR failed to converge after %u iterations: branch %3d - %3d "
        "(old: %.12f, new: %.12f)\n",
        params->max_newton_iters, tr_p->clv_index, tr_p->back->clv_index,
        tr_p->length, xres);

    if (!check_loglh_improvement(params->opt_method))
      xres = tr_p->length;

    pllmod_reset_error();
  }

  if (pll_errno)
    return PLL_FAILURE;

  assert(xres >= xmin && xres <= xmax);

  if (fabs(tr_p->length - xres) > 1e-10)
  {
    /* update branch length in the tree structure */
    tr_p->length = tr_p->back->length = xres;

    /* update pmatrix for the new branch length */
    if (keep_update)
    {
      update_prob_matrices(params->partitions, params->partition_count,
                           params->params_indices, params->brlen_scalers, tr_p);
    }

    if (check_loglh_improvement(params->opt_method))
    {
      assert(keep_update);

      /* check and compare likelihood */
      double eval_loglikelihood = pllmod_opt_compute_edge_loglikelihood_multi(
                                                        params->partitions,
                                                        params->partition_count,
                                                        tr_p->clv_index,
                                                        tr_p->scaler_index,
                                                        tr_p->back->clv_index,
                                                        tr_p->back->scaler_index,
                                                        tr_p->pmatrix_index,
                                                        params->params_indices,
                                                        NULL,
                                                        params->parallel_context,
                                                        params->parallel_reduce_cb);

      /* check if the optimal found value improves the likelihood score */
      if (eval_loglikelihood >= *loglikelihood_score)
      {
        DBG("ACCEPT: new BL: %.12lf, old BL: %.12lf, new LH: %.9lf, old LH: %.9lf\n",
               xres, xorig, eval_loglikelihood, *loglikelihood_score);

        /* fix new score */
        *loglikelihood_score = eval_loglikelihood;
      }
      else
      {
        DBG("REVERT: new BL: %.12lf, old BL: %.12lf, new LH: %.9lf, old LH: %.9lf\n",
               xres, xorig, eval_loglikelihood, *loglikelihood_score);

        tr_p->length = tr_p->back->length = xorig;

        /* reset branch length */
        update_prob_matrices(params->partitions, params->partition_count,
                             params->params_indices, params->brlen_scalers, tr_p);
      }
    }
  }

#ifdef DEBUG
  {
    DBG(" Optimized branch %3d - %3d (%.12f -> %.12f)\n",
        tr_p->clv_index, tr_p->back->clv_index, xorig, tr_p->length);

    double new_loglh = pllmod_opt_compute_edge_loglikelihood_multi(
                                                      params->partitions,
                                                      params->partition_count,
                                                      tr_p->clv_index,
                                                      tr_p->scaler_index,
                                                      tr_p->back->clv_index,
                                                      tr_p->back->scaler_index,
                                                      tr_p->pmatrix_index,
                                                      params->params_indices,
                                                      NULL,
                                                      params->parallel_context,
                                                      params->parallel_reduce_cb);

    DBG(" New loglH: %.12f\n", new_loglh);
  }
#endif

  /* update children */
  if (radius && tr_q && tr_z)
  {
    /* update children 'Q'
     * CLV at P is recomputed with children P->back and Z->back
     * Scaler is updated by subtracting Q->back and adding P->back
     */
    update_partials_and_scalers(params->partitions,
                                params->partition_count,
                                tr_q,
                                tr_p,
                                tr_z);

    /* eval */
    pll_newton_tree_params_multi_t params_cpy;
    memcpy(&params_cpy, params, sizeof(pll_newton_tree_params_multi_t));
    params_cpy.tree = tr_q->back;
    if (!recomp_iterative_multi (&params_cpy,
                           radius-1,
                           loglikelihood_score,
                           keep_update))
      return PLL_FAILURE;

    /* update children 'Z'
     * CLV at P is recomputed with children P->back and Q->back
     * Scaler is updated by subtracting Z->back and adding Q->back
     */
    update_partials_and_scalers(params->partitions,
                                params->partition_count,
                                tr_z,
                                tr_q,
                                tr_p);

   /* eval */
    params_cpy.tree = tr_z->back;
    if (!recomp_iterative_multi (&params_cpy,
                           radius-1,
                           loglikelihood_score,
                           keep_update))
      return PLL_FAILURE;

    /* reset to initial state
     * CLV at P is recomputed with children Q->back and Z->back
     * Scaler is updated by subtracting P->back and adding Z->back
     */
    update_partials_and_scalers(params->partitions,
                                params->partition_count,
                                tr_p,
                                tr_z,
                                tr_q);
  }

  return PLL_SUCCESS;

} /* recomp_iterative */

/**
 * Optimize branch lengths locally around a given edge using Newton-Raphson
 * minimization algorithm on a multiple partition.
 *
 * Check `pllmod_opt_optimize_branch_lengths_local` documentation.
 *
 * @param[in,out]  partitions list of partitions
 * @param  partition_count    number of partitions in `partitions`
 * @param[in,out]  tree       the PLL unrooted tree structure
 * @param  params_indices     the indices of the parameter sets
 * @param  precomp_buffers    buffer for sumtable (NULL=allocate internally)
 * @param  brlen_buffers      buffer for branch lengths (NULL=allocate internally)
 * @param  brlen_scalers      branch length scalers
 * @param  branch_length_min  lower bound for branch lengths
 * @param  branch_length_max  upper bound for branch lengths
 * @param  tolerance          tolerance for Newton-Raphson algorithm
 * @param  smoothings         number of iterations over the branches
 * @param  radius             radius from the virtual root
 * @param  keep_update        if true, branch lengths are iteratively updated in the tree structure
 * @param  opt_method         optimization method to use (see PLLMOD_OPT_BLO_* constants)
 * @param  parallel_context   context for parallel computation
 * @param  parallel_reduce_cb callback function for parallel reduction
 *
 * @return                   the likelihood score after optimizing branch lengths
 */
PLL_EXPORT double pllmod_opt_optimize_branch_lengths_local_multi (
                                              pll_partition_t ** partitions,
                                              size_t partition_count,
                                              pll_unode_t * tree,
                                              unsigned int ** params_indices,
                                              double ** precomp_buffers,
                                              double ** brlen_buffers,
                                              double * brlen_scalers,
                                              double branch_length_min,
                                              double branch_length_max,
                                              double tolerance,
                                              int smoothings,
                                              int radius,
                                              int keep_update,
                                              int opt_method,
                                              int brlen_linkage,
                                              void * parallel_context,
                                              void (*parallel_reduce_cb)(void *,
                                                                         double *,
                                                                         size_t,
                                                                         int))
{
  unsigned int iters;
  double loglikelihood = 0.0, new_loglikelihood;
  size_t p;
  double result = (double) PLL_FAILURE;

  pllmod_reset_error();

  /**
   * preconditions:
   *    (1) CLVs must be updated towards 'tree'
   *    (2) Pmatrix indices must be **unique** for each branch
   */

  if (opt_method == PLLMOD_OPT_BLO_NEWTON_FALLBACK ||
      opt_method == PLLMOD_OPT_BLO_NEWTON_GLOBAL)
  {
    pllmod_set_error(PLLMOD_ERROR_NOT_IMPLEMENTED,
                     "Optimization method not implemented: "
                     "NEWTON_FALLBACK, NEWTON_GLOBAL");
    return (double)PLL_FAILURE;
  }


  if (brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED)
  {
    pllmod_set_error(PLLMOD_ERROR_NOT_IMPLEMENTED,
                     "Optimization of unlinked branch length is not implemented yet.");
    return (double)PLL_FAILURE;
  }

  if (radius < PLLMOD_OPT_BRLEN_OPTIMIZE_ALL)
  {
    pllmod_set_error(PLLMOD_OPT_ERROR_NEWTON_BAD_RADIUS,
                     "Invalid radius for branch length optimization");
    return (double)PLL_FAILURE;
  }

  /* get the initial likelihood score */
  loglikelihood = pllmod_opt_compute_edge_loglikelihood_multi (
                                                      partitions,
                                                      partition_count,
                                                      tree->back->clv_index,
                                                      tree->back->scaler_index,
                                                      tree->clv_index,
                                                      tree->scaler_index,
                                                      tree->pmatrix_index,
                                                      params_indices,
                                                      NULL,
                                                      parallel_context,
                                                      parallel_reduce_cb);

  DBG("\nStarting BLO_multi: radius: %u, max_iters: %u, old LH: %.9f\n",
      radius, smoothings, loglikelihood);

  /* set parameters for N-R optimization */
  pll_newton_tree_params_multi_t params;
  params.partitions        = partitions;
  params.partition_count   = partition_count;
  params.tree              = tree;
  params.params_indices    = params_indices;
  params.branch_length_min = (branch_length_min>0)?
                              branch_length_min:
                              PLLMOD_OPT_MIN_BRANCH_LEN;
  params.branch_length_max = (branch_length_max>0)?
                              branch_length_max:
                              PLLMOD_OPT_MAX_BRANCH_LEN;
  params.tolerance         = (branch_length_min>0)?
                              branch_length_min/10.0:
                              PLLMOD_OPT_TOL_BRANCH_LEN;
  params.precomp_buffers   = precomp_buffers;
  params.brlen_buffers     = brlen_buffers;
  params.brlen_scalers     = brlen_scalers;
  params.opt_method        = opt_method;
  params.brlen_linkage     = brlen_linkage;
  params.max_newton_iters  = 30;

  params.parallel_context = parallel_context;
  params.parallel_reduce_cb = parallel_reduce_cb;

  /* allocate the sumtable if needed */
  if (!allocate_buffers(&params))
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for bl opt variables");
    goto cleanup;
  }

  iters = (unsigned int) smoothings;
  while (iters)
  {
    new_loglikelihood = loglikelihood;

    /* iterate on first edge */
    params.tree = tree;
    if (!recomp_iterative_multi (&params, radius, &new_loglikelihood, keep_update))
      goto cleanup;

    if (radius)
    {
      /* iterate on second edge */
      params.tree = tree->back;
      if (!recomp_iterative_multi (&params, radius-1, &new_loglikelihood, keep_update))
        goto cleanup;
    }

    /* compute likelihood after optimization */
    new_loglikelihood = pllmod_opt_compute_edge_loglikelihood_multi (
                                                        partitions,
                                                        partition_count,
                                                        tree->back->clv_index,
                                                        tree->back->scaler_index,
                                                        tree->clv_index,
                                                        tree->scaler_index,
                                                        tree->pmatrix_index,
                                                        params_indices,
                                                        NULL,
                                                        parallel_context,
                                                        parallel_reduce_cb);

    DBG("BLO_multi: iteration %u, old LH: %.9f, new LH: %.9f\n",
        (unsigned int) smoothings - iters, loglikelihood, new_loglikelihood);

    if (new_loglikelihood - loglikelihood > new_loglikelihood * BETTER_LL_TRESHOLD)
    {
      iters --;

      /* check convergence */
      if (fabs (new_loglikelihood - loglikelihood) < tolerance) iters = 0;

      loglikelihood = new_loglikelihood;
    }
    else
    {
      if (params.opt_method == PLLMOD_OPT_BLO_NEWTON_SAFE ||
          params.opt_method == PLLMOD_OPT_BLO_NEWTON_OLDSAFE)
        assert(new_loglikelihood - loglikelihood > new_loglikelihood * BETTER_LL_TRESHOLD);
      else if (opt_method == PLLMOD_OPT_BLO_NEWTON_FALLBACK)
      {
        // reset branch lenghts
        params.opt_method = PLLMOD_OPT_BLO_NEWTON_SAFE;
        iters = smoothings;
      }
      else
      {
        pllmod_set_error(PLLMOD_OPT_ERROR_NEWTON_WORSE_LK,
                         "BL opt converged to a worse likelihood score by %.15f units",
                         new_loglikelihood - loglikelihood);
        goto cleanup;
      }
    }

  }  // while

  result = -1*loglikelihood;

cleanup:
  /* deallocate sumtable */
  if (!precomp_buffers)
  {
    for (p = 0; p < partition_count; ++p)
    {
      if (params.precomp_buffers[p])
        free(params.precomp_buffers[p]);
    }
    pll_aligned_free(params.precomp_buffers);
  }

  if (params.brlen_buffers && !brlen_buffers)
  {
    if (params.brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED)
    {
      for (p = 0; p < partition_count; ++p)
      {
        if (params.brlen_buffers[p])
          free(params.brlen_buffers[p]);
      }
    }
    else
      free(params.brlen_buffers[0]);
    pll_aligned_free(params.brlen_buffers);
  }

  return result;
} /* pllmod_opt_optimize_branch_lengths_local */
