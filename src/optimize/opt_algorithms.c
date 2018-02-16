/*
 Copyright (C) 2015 Diego Darriba

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
#include "pll_optimize.h"
#include "lbfgsb/lbfgsb.h"
#include "../pllmod_common.h"

/**
 * @file opt_algorithms.c
 *
 * @brief Core optimization algorithms
 *
 * This file implements low-level optimization algorithms.
 * Algorithms here are not specific for phylogenetic analysis and can be used
 * given the proper parameters. In general, boundaries for the variable/s and
 * a target function that returns a score (to be minimized) given parameter
 * value/s.
 *
 * @author Diego Darriba
 * @author Alexey Kozlov
 */

static inline int is_nan(double v)
{
  return v!=v;
}

static inline int d_equals(double a, double b)
{
  return (fabs(a-b) < 1e-10);
}

/******************************************************************************/
/* NEWTON-RAPHSON OPTIMIZATION */
/******************************************************************************/

/**
 * Minimize a function using Newton-Raphson.
 *
 * The target function must compute the derivatives at a certain point,
 * and it requires 4 parameters: (1) custom data (if needed), (2) the value at
 * which derivatives are computed, and (3,4) lower and upper bounds.
 *
 * @param  x1         lower bound
 * @param  xguess     first guess for the free variable
 * @param  x2         upper bound
 * @param  tolerance  tolerance of the minimization method
 * @param  max_iters  maximum number of iterations (bounds the effect of slow convergence)
 * @param  params     custom parameters required by the target function
 * @param  deriv_func target function
 *
 * @return            the parameter value that minimizes the function in [x1,x2]
 */
PLL_EXPORT double pllmod_opt_minimize_newton(double xmin,
                                             double xguess,
                                             double xmax,
                                             double tolerance,
                                             unsigned int max_iters,
                                             void * params,
                                             void (deriv_func)(void *,
                                                          double,
                                                          double *, double *))
{
  unsigned int i;
  double dx, f, df;
  double xh, xl, x;

  double dxmax = xmax / max_iters;

  /* reset errno */
  pll_errno = 0;

  x = PLL_MAX(PLL_MIN(xguess, xmax), xmin);

  xl = xmin;
  xh = xmax;

  for (i = 1; i <= max_iters; i++)
  {
    deriv_func((void *)params, x, &f, &df);

    if (!isfinite(f) || !isfinite(df))
    {
      DBG("[%d][NR deriv] BL=%.9f   f=%.12f  df=%.12f\n",
          i, x, f, df);
      pllmod_set_error(PLLMOD_OPT_ERROR_NEWTON_DERIV,
                       "Wrong likelihood derivatives");
      return PLL_FAILURE;
    }

    if (df > 0.0)
    {
      if (fabs (f) < tolerance)
        return x;

      if (f < 0.0)
        xl = x;
      else
        xh = x;

      dx = -1 * f / df;
    }
    else
    {
//          dx = xl[root_num] + 0.5 * (xh[root_num] - xl[root_num]) - x[root_num];
      dx = -1 * f / fabs(df);
    }

    dx = PLL_MAX(PLL_MIN(dx, dxmax), -dxmax);

    if (x+dx < xl) dx = xl-x;
    if (x+dx > xh) dx = xh-x;

//    if ((x+dx < xl) || (x+dx > xh))
//    {
//      // reset to the middle of the current interval
//      dx = xl + 0.5 * (xh - xl) - x;
//    }

    if (fabs (dx) < tolerance)
      return x;

    DBG("[%d][NR deriv] BL=%.9f   f=%.12f  df=%.12f  nextBL=%.9f\n",
        i, x, f, df, x+dx);

    x += dx;

    x = PLL_MAX(PLL_MIN(x, xmax), xmin);
  }

  pllmod_set_error(PLLMOD_OPT_ERROR_NEWTON_LIMIT,
                   "Exceeded maximum number of iterations");
  return PLL_FAILURE;
}


/**
 * Minimize a function using Newton-Raphson - LEGACY version from IQTree
 *
 * The target function must compute the derivatives at a certain point,
 * and it requires 4 parameters: (1) custom data (if needed), (2) the value at
 * which derivatives are computed, and (3,4) lower and upper bounds.
 *
 * @param  x1         lower bound
 * @param  xguess     first guess for the free variable
 * @param  x2         upper bound
 * @param  tolerance  tolerance of the minimization method
 * @param  max_iters  maximum number of iterations (bounds the effect of slow convergence)
 * @param  params     custom parameters required by the target function
 * @param  deriv_func target function
 *
 * @return            the parameter value that minimizes the function in [x1,x2]
 */
PLL_EXPORT double pllmod_opt_minimize_newton_old(double x1,
                                                  double xguess,
                                                  double x2,
                                                  double tolerance,
                                                  unsigned int max_iters,
                                                  void * params,
                                                  void (deriv_func)(void *,
                                                                    double,
                                                                    double *, double *))
{
  unsigned int i;
  double df, dx, dxold, f;
  double temp, xh, xl, rts, rts_old = 0.0;

  /* reset errno */
  pll_errno = 0;

  rts = xguess;
  if (rts < x1)
    rts = x1;
  if (rts > x2)
    rts = x2;

  deriv_func((void *)params, rts, &f, &df);

  DBG("[NR deriv] BL=%.9f   f=%.12f  df=%.12f  nextBL=%.9f\n", rts, f, df, rts-f/fabs(df));
  if (!isfinite(f) || !isfinite(df))
  {
    pllmod_set_error(PLLMOD_OPT_ERROR_NEWTON_DERIV,
                     "wrong likelihood derivatives");
    return (double) -INFINITY;
  }
  if (df >= 0.0 && fabs (f) < tolerance)
    return rts;
  if (f < 0.0)
  {
    xl = rts;
    xh = x2;
  }
  else
  {
    xh = rts;
    xl = x1;
  }

  dx = dxold = fabs (xh - xl);
  for (i = 1; i <= max_iters; i++)
  {
    rts_old = rts;

    if ((df <= 0.0) // function is concave
    || (((rts - xh) * df - f) * ((rts - xl) * df - f) >= 0.0) // out of bound
        )
    {
      dxold = dx;
      dx = 0.5 * (xh - xl);
      rts = xl + dx;
      if (xl == rts)
        return rts;
    }
    else
    {
      dxold = dx;
      dx = f / df;
      temp = rts;
      rts -= dx;
      if (temp == rts)
        return rts;
    }

    if (fabs (dx) < tolerance)
      return rts_old;

    if (i == max_iters)
      break;

    if (rts < x1) rts = x1;

    deriv_func((void *)params, rts, &f, &df);

    DBG("[%d][NR deriv] BL=%.9f   f=%.12f  df=%.12f  nextBL=%.9f\n", i, rts, f, df, rts-f/df);

    if (!isfinite(f) || !isfinite(df))
    {
      pllmod_set_error(PLLMOD_OPT_ERROR_NEWTON_DERIV,
                       "Wrong likelihood derivatives [it]");
      return (double) -INFINITY;
    }

    if (df > 0.0 && fabs (f) < tolerance)
    {
      return rts;
    }

    if (f < 0.0)
      xl = rts;
    else
      xh = rts;
  }

  pllmod_set_error(PLLMOD_OPT_ERROR_NEWTON_LIMIT,
                   "Exceeded maximum number of iterations");
  return rts_old;
}


/******************************************************************************/
/* L-BFGS-B OPTIMIZATION */
/******************************************************************************/

/**
 * Minimize a multi-parameter function using L-BFGS-B.
 *
 * The target function must compute the score at a certain state,
 * and it requires 2 parameters: (1) custom data (if needed),
 * and (2) the values at which score is computed.
 *
 * @param  x[in,out]   first guess and result of the minimization process
 * @param  xmin        lower bound for each of the variables
 * @param  xmax        upper bound for each of the variables
 * @param  bound       bound type (PLL_LBFGSB_BOUND_[NONE|LOWER|UPPER|BOTH]
 * @param  n           number of variables
 * @param  factr       convergence tolerance for L-BFGS-B relative to machine epsilon
 * @param  pgtol       absolute gradient tolerance for L-BFGS-B
 * @param  params      custom parameters required by the target function
 * @param  target_funk target function
 *
 * `factr` is a double precision variable. The iteration will stop when
 * (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
 * where epsmch is the machine epsilon
 *
 * `pgtol` is a double precision variable. The iteration will stop when
 * max{|proj g_i | i = 1, ..., n} <= pgtol
 * where pg_i is the ith component of the projected gradient.
 *
 * @return             the minimal score found
 */
PLL_EXPORT double pllmod_opt_minimize_lbfgsb (double * x,
                                             double * xmin,
                                             double * xmax,
                                             int * bound,
                                             unsigned int n,
                                             double factr,
                                             double pgtol,
                                             void * params,
                                             double (*target_funk)(
                                                     void *,
                                                     double *))
{
  unsigned int i;

  /* L-BFGS-B parameters */
  //  double initial_score;
  int max_corrections;
  double score = 0;
  double *g, *wa;
  int *iwa;

  int taskValue;
  int *task = &taskValue;

  int csaveValue;
  int *csave = &csaveValue;
  double dsave[29];
  int isave[44];
  logical lsave[4];

  int iprint = -1;

  max_corrections = 5;

  /* reset errno */
  pll_errno = 0;

  g = (double *) calloc ((size_t) n, sizeof(double));

  /*     We start the iteration by initializing task. */
  *task = (int) START;

  iwa = (int *) calloc (3 * (size_t) n, sizeof(int));
  wa = (double *) calloc (
      (2 * (size_t)max_corrections + 5) * (size_t)n
          + 12 * (size_t)max_corrections * ((size_t)max_corrections + 1),
      sizeof(double));

  if (!(wa && iwa && g))
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for l-bfgs-b variables");
    if (g)
      free (g);
    if (iwa)
      free (iwa);
    if (wa)
      free (wa);
    return (double) -INFINITY;
  }

//  double initial_score = target_funk (params, x);
  int continue_opt = 1;
  while (continue_opt)
  {
    /*     This is the call to the L-BFGS-B code. */
    setulb ((int *)&n, &max_corrections, x, xmin, xmax,
            bound, &score, g, &factr, &pgtol, wa, iwa,
            task, &iprint, csave, lsave, isave, dsave);
    if (IS_FG(*task))
    {
      /*
       * the minimization routine has returned to request the
       * function f and gradient g values at the current x.
       * Compute function value f for the sample problem.
       */

      score = target_funk(params, x);

      if (is_nan(score) || d_equals(score, (double) -INFINITY))
        break;

      double h, temp;
      for (i = 0; i < n; i++)
      {
        temp = x[i];
        h = PLL_LBFGSB_ERROR * fabs (temp);
        if (h < 1e-12)
          h = PLL_LBFGSB_ERROR;

        x[i] = temp + h;
        h = x[i] - temp;
        double lnderiv = target_funk(params, x);

        g[i] = (lnderiv - score) / h;

        /* reset variable */
        x[i] = temp;
      }
    }
    else if (*task != NEW_X)
      continue_opt = 0;
  }

  /* fix optimal parameters */
  score = target_funk(params, x);

  free (iwa);
  free (wa);
  free (g);

  if (is_nan(score))
  {
    score = (double) -INFINITY;
    /* set errno only if it was not set by some inner function */
    if (!pll_errno)
    {
      pllmod_set_error(PLLMOD_OPT_ERROR_LBFGSB_UNKNOWN, "Unknown LBFGSB error");
    }
  }

  return score;
} /* pllmod_opt_minimize_lbfgsb */

struct bfgs_multi_opt
{
  unsigned int n;
  double * x;
  double * xmin;
  double * xmax;
  int * bound;
  double factr;
  double pgtol;

  int iprint;
  int max_corrections;
  int task;

  double *g, *wa;
  int *iwa;
  int csave;
  double dsave[29];
  int isave[44];
  logical lsave[4];

  double score;
  double h, temp;
};

static int init_bfgs_opt(struct bfgs_multi_opt * opt, unsigned int n, double * x,
                         double * xmin, double * xmax, int * bound, double factr,
                         double pgtol)
{
  /* We start the iteration by initializing task. */
  opt->task = (int) START;
  opt->score = 0;

  // some magic numbers
  opt->iprint = -1;
  opt->max_corrections = 5;

  opt->n = n;
  opt->x = x;
  opt->xmin = xmin;
  opt->xmax = xmax;
  opt->bound = bound;
  opt->factr = factr;
  opt->pgtol = pgtol;

  /* allocate memory */
  opt->g = (double *) calloc ((size_t) n, sizeof(double));

  opt->iwa = (int *) calloc (3 * (size_t) n, sizeof(int));
  opt->wa = (double *) calloc (
      (2 * (size_t)opt->max_corrections + 5) * (size_t)n
          + 12 * (size_t)opt->max_corrections * ((size_t)opt->max_corrections + 1),
      sizeof(double));

  if (!opt->g || !opt->iwa || !opt->wa)
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for l-bfgs-b variables");
    return PLL_FAILURE;
  }

  return PLL_SUCCESS;
}

static void destroy_bfgs_opt(struct bfgs_multi_opt * opt)
{
  if (opt)
  {
    if (opt->g)
      free(opt->g);
    if (opt->iwa)
      free(opt->iwa);
    if (opt->wa)
      free(opt->wa);
    free(opt);
  }
}

static int setulb_multi(struct bfgs_multi_opt * opt)
{
  return setulb ((int *)&opt->n, &opt->max_corrections, opt->x, opt->xmin,
                 opt->xmax, opt->bound, &opt->score, opt->g, &opt->factr,
                 &opt->pgtol, opt->wa, opt->iwa, &opt->task, &opt->iprint,
                 &opt->csave, opt->lsave, opt->isave, opt->dsave);
}

PLL_EXPORT double pllmod_opt_minimize_lbfgsb_multi(unsigned int xnum,
                                                   double ** x,
                                                   double ** xmin,
                                                   double ** xmax,
                                                   int ** bound,
                                                   unsigned int * n,
                                                   unsigned int nmax,
                                                   double factr,
                                                   double pgtol,
                                                   void * params,
                                                   double (*target_funk)(void *,
                                                                         double **,
                                                                         double *,
                                                                         int *))
{
  unsigned int i, p;

  double score = (double) -INFINITY;

  double * lh_old = (double *) calloc ((size_t) xnum, sizeof(double));
  double * lh_new = (double *) calloc ((size_t) xnum, sizeof(double));
  int * converged = (int *) calloc ((size_t) xnum+1, sizeof(int));
  int * skip = (int *) calloc ((size_t) xnum+1, sizeof(int));

  struct bfgs_multi_opt ** opts = (struct bfgs_multi_opt **)
                           calloc((size_t) xnum, sizeof(struct bfgs_multi_opt *));

  if (!lh_old || !lh_new || !converged || !skip || !opts)
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for l-bfgs-b variables");
    goto cleanup;
  }

  for (p = 0; p < xnum; p++)
  {
    if(!x[p])
    {
      /* this is a remote partition - forget about it */
      converged[p] = skip[p] = 1;
      continue;
    }

    opts[p] = (struct bfgs_multi_opt *) calloc(1, sizeof(struct bfgs_multi_opt));
    if (!opts[p])
    {
      pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                       "Cannot allocate memory for l-bfgs-b variables");
      goto cleanup;
    }

    if (!init_bfgs_opt(opts[p], n[p], x[p], xmin[p], xmax[p], bound[p], factr,
                       pgtol))
      goto cleanup;
  }

  /* reset errno */
  pll_errno = 0;

  int continue_opt = 1;
  while (continue_opt)
  {
    continue_opt = 0;
    int all_skip = 1;
    for (p = 0; p < xnum; p++)
    {
      if (converged[p])
      {
        skip[p] = 1;
        continue;
      }

      assert(opts[p]);

      /*     This is the call to the L-BFGS-B code. */
      setulb_multi(opts[p]);

      skip[p] = !IS_FG(opts[p]->task);
      converged[p] = skip[p] && (opts[p]->task != NEW_X);
      all_skip &= skip[p];
    }

    /* check if ALL partitions have converged, including remote ones */
    target_funk (params, NULL, NULL, converged);
    continue_opt = !converged[xnum];
    target_funk (params, NULL, NULL, skip);
    all_skip = converged[xnum];

    if (!all_skip && continue_opt)
    {
      /*
       * the minimization routine has returned to request the
       * function f and gradient g values at the current x.
       * Compute function value f for the sample problem.
       */
      score = target_funk (params, x, lh_old, skip);

      if (is_nan(score) || d_equals(score, (double) -INFINITY))
        break;

      for (p = 0; p < xnum; p++)
      {
        if (!skip[p])
          opts[p]->score = lh_old[p];
      }

      for (i = 0; i < nmax; i++)
      {
        for (p = 0; p < xnum; p++)
        {
          if (i >= n[p])
          {
            /* this partition has less parameters than max -> skip it */
            skip[p] = 1;
          }

          if (skip[p])
            continue;

          opts[p]->temp = x[p][i];
          opts[p]->h = PLL_LBFGSB_ERROR * fabs (opts[p]->temp);
          if (opts[p]->h < 1e-12)
            opts[p]->h = PLL_LBFGSB_ERROR;

          x[p][i] = opts[p]->temp + opts[p]->h;
          opts[p]->h = x[p][i] - opts[p]->temp;
        }

        score = target_funk(params, x, lh_new, skip);

        for (p = 0; p < xnum; p++)
        {
          if (skip[p])
            continue;

          assert(opts[p]);

          /* compute partial derivative */
          opts[p]->g[i] = (lh_new[p] - lh_old[p]) / opts[p]->h;

          /* reset variable */
          x[p][i] = opts[p]->temp;

          DBG("P%d  lh_old: %f, lh_new: %f\n", p, lh_old[p], lh_new[p]);
        }
      }
    }
  }

  /* fix optimal parameters */
  score = target_funk (params, x, NULL, NULL);

cleanup:
  if (lh_old)
    free(lh_old);
  if (lh_new)
    free(lh_new);
  if (converged)
    free(converged);
  if (skip)
    free(skip);
  if (opts)
  {
    for (p = 0; p < xnum; p++)
      destroy_bfgs_opt(opts[p]);
    free(opts);
  }

  if (is_nan(score))
  {
    score = (double) -INFINITY;
    /* set errno only if it was not set by some inner function */
    if (!pll_errno)
    {
      pllmod_set_error(PLLMOD_OPT_ERROR_LBFGSB_UNKNOWN, "Unknown LBFGSB error");
    }
  }

  return score;
} /* pllmod_opt_minimize_lbfgsb */

/******************************************************************************/
/* BRENT'S OPTIMIZATION */
/******************************************************************************/

#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-7
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

struct opt_params
{
  double startx;
  double fstartx;

  double tol;
  double *foptx;
  double *f2optx;

  int iter;
  double a;
  double b;
  double d;
  double etemp;
  double fu;
  double fv;
  double fw;
  double fx;
  double p;
  double q;
  double r;
  double tol1;
  double tol2;
  double u;
  double v;
  double w;
  double x;
  double xm;
  double xw;
  double wv;
  double vx;
  double e;
};

struct brent_wrapper_params
{
  double (*target_funk)(void *, double);
  void * params;
};

static int brent_opt_init (double ax, double bx, double cx, double tol,
                         double *foptx, double *f2optx, double fax,
                         double fbx, double fcx,
                         struct opt_params * bp)
{
  memset(bp, 0, sizeof(struct opt_params));
  double etemp;

  bp->tol = tol;
  bp->foptx = foptx;
  bp->f2optx = f2optx;

  bp->a = (ax < cx ? ax : cx);
  bp->b = (ax > cx ? ax : cx);
  bp->startx = bp->x = bx;
  bp->fstartx = bp->fx = fbx;
  if (fax < fcx)
  {
    bp->w = ax;
    bp->fw = fax;
    bp->v = cx;
    bp->fv = fcx;
  }
  else
  {
    bp->w = cx;
    bp->fw = fcx;
    bp->v = ax;
    bp->fv = fax;
  }

  /* pre-loop iteration 0 */
  bp->iter = 1;

  bp->xm = 0.5 * (bp->a + bp->b);
  bp->tol2 = 2.0 * (bp->tol1 = bp->tol * fabs (bp->x) + ZEPS);
  if (fabs (bp->x - bp->xm) <= (bp->tol2 - 0.5 * (bp->b - bp->a)))
  {
    if (bp->foptx)
      *bp->foptx = bp->fx;
    bp->xw = bp->x - bp->w;
    bp->wv = bp->w - bp->v;
    bp->vx = bp->v - bp->x;
    if (bp->f2optx)
    {
      *bp->f2optx = 2.0 * (bp->fv * bp->xw + bp->fx * bp->wv + bp->fw * bp->vx)
          / (bp->v * bp->v * bp->xw + bp->x * bp->x * bp->wv + bp->w * bp->w * bp->vx);
    }
    return PLL_FAILURE;
  }

  if (fabs (bp->e) > bp->tol1)
  {
    bp->r = (bp->x - bp->w) * (bp->fx - bp->fv);
    bp->q = (bp->x - bp->v) * (bp->fx - bp->fw);
    bp->p = (bp->x - bp->v) * bp->q - (bp->x - bp->w) * bp->r;
    bp->q = 2.0 * (bp->q - bp->r);
    if (bp->q > 0.0)
      bp->p = -bp->p;
    bp->q = fabs (bp->q);
    etemp = bp->e;
    bp->e = bp->d;
    if (fabs (bp->p) >= fabs (0.5 * bp->q * etemp) || bp->p <= bp->q * (bp->a - bp->x)
        || bp->p >= bp->q * (bp->b - bp->x))
      bp->d = CGOLD * (bp->e = (bp->x >= bp->xm ? bp->a - bp->x : bp->b - bp->x));
    else /*TODO:CONTINUE HERE!!!!!!!!!!!!!!!!!!!!!!!!!*/
    {
      bp->d = bp->p / bp->q;
      bp->u = bp->x + bp->d;
      if (bp->u - bp->a < bp->tol2 || bp->b - bp->u < bp->tol2)
        bp->d = SIGN(bp->tol1, bp->xm - bp->x);
    }
  }
  else
  {
    bp->d = CGOLD * (bp->e = (bp->x >= bp->xm ? bp->a - bp->x : bp->b - bp->x));
  }

  bp->u = (fabs (bp->d) >= bp->tol1 ? bp->x + bp->d : bp->x + SIGN(bp->tol1, bp->d));
  return PLL_SUCCESS;
}

static int brent_opt_post_loop (struct opt_params * bp)
{
  double etemp;

  /* post-loop iteration i */

  if (bp->fu <= bp->fx)
  {
    if (bp->u >= bp->x)
      bp->a = bp->x;
    else
      bp->b = bp->x;

    SHFT(bp->v, bp->w, bp->x, bp->u)
    SHFT(bp->fv, bp->fw, bp->fx, bp->fu)
  }
  else
  {
    if (bp->u < bp->x)
      bp->a = bp->u;
    else
      bp->b = bp->u;
    if (bp->fu <= bp->fw || bp->w == bp->x)
    {
      bp->v = bp->w;
      bp->w = bp->u;
      bp->fv = bp->fw;
      bp->fw = bp->fu;
    }
    else if (bp->fu <= bp->fv || bp->v == bp->x || bp->v == bp->w)
    {
      bp->v = bp->u;
      bp->fv = bp->fu;
    }
  }

  /* pre-loop iteration i+1*/

  ++bp->iter;

  bp->xm = 0.5 * (bp->a + bp->b);
  bp->tol2 = 2.0 * (bp->tol1 = bp->tol * fabs (bp->x) + ZEPS);
  if (fabs (bp->x - bp->xm) <= (bp->tol2 - 0.5 * (bp->b - bp->a)))
  {
    if (bp->foptx)
      *bp->foptx = bp->fx;
    bp->xw = bp->x - bp->w;
    bp->wv = bp->w - bp->v;
    bp->vx = bp->v - bp->x;
    if (bp->f2optx)
    {
      *bp->f2optx = 2.0 * (bp->fv * bp->xw + bp->fx * bp->wv + bp->fw * bp->vx)
          / (bp->v * bp->v * bp->xw + bp->x * bp->x * bp->wv + bp->w * bp->w * bp->vx);
    }
    return PLL_FAILURE;
  }

  if (fabs (bp->e) > bp->tol1)
  {
    bp->r = (bp->x - bp->w) * (bp->fx - bp->fv);
    bp->q = (bp->x - bp->v) * (bp->fx - bp->fw);
    bp->p = (bp->x - bp->v) * bp->q - (bp->x - bp->w) * bp->r;
    bp->q = 2.0 * (bp->q - bp->r);
    if (bp->q > 0.0)
      bp->p = -bp->p;
    bp->q = fabs (bp->q);
    etemp = bp->e;
    bp->e = bp->d;
    if (fabs (bp->p) >= fabs (0.5 * bp->q * etemp) || bp->p <= bp->q * (bp->a - bp->x)
        || bp->p >= bp->q * (bp->b - bp->x))
      bp->d = CGOLD * (bp->e = (bp->x >= bp->xm ? bp->a - bp->x : bp->b - bp->x));
    else
    {
      bp->d = bp->p / bp->q;
      bp->u = bp->x + bp->d;
      if (bp->u - bp->a < bp->tol2 || bp->b - bp->u < bp->tol2)
        bp->d = SIGN(bp->tol1, bp->xm - bp->x);
    }
  }
  else
  {
    bp->d = CGOLD * (bp->e = (bp->x >= bp->xm ? bp->a - bp->x : bp->b - bp->x));
  }

  bp->u = (fabs (bp->d) >= bp->tol1 ? bp->x + bp->d : bp->x + SIGN(bp->tol1, bp->d));

  return PLL_SUCCESS;
}

double target_funk_wrapper(void * params,
                           double * xopt,
                           double * fxopt,
                           int * converged)
{
  struct brent_wrapper_params * wrap_params =
      (struct brent_wrapper_params *) params;

  *fxopt = wrap_params->target_funk(wrap_params->params, *xopt);

  return *fxopt;
}

static int brent_opt_alt (int xnum,
                          int * opt_mask,
                          double * xmin,
                          double * xguess,
                          double * xmax,
                          double xtol,
                          double * xopt,
                          double * fx,
                          double * f2x,
                          void * params,
                          double (*target_funk)(
                                  void *,
                                  double *,
                                  double *,
                                  int *),
                          int global_range)
{
  struct opt_params * brent_params = (struct opt_params *)
      calloc(xnum, sizeof(struct opt_params));

  double * ax = (double *) calloc(xnum, sizeof(double));
  double * cx = (double *) calloc(xnum, sizeof(double));
  double * fa = (double *) calloc(xnum, sizeof(double));
  double * fb = (double *) calloc(xnum, sizeof(double));
  double * fc = (double *) calloc(xnum, sizeof(double));
  double * fxmin = (double *) calloc(xnum, sizeof(double));
  double * fxmax = (double *) calloc(xnum, sizeof(double));
  int * init_ok = (int *) calloc(xnum+1, sizeof(int));


  double * l_xmin = NULL;
  double * l_xmax = NULL;

  int i;
  int iterate = 1;

  if (global_range)
  {
    l_xmin = (double *) calloc(xnum, sizeof(double));
    l_xmax = (double *) calloc(xnum, sizeof(double));
    for (i = 0; i < xnum; ++i)
    {
      l_xmin[i] = *xmin;
      l_xmax[i] = *xmax;
    }
  }
  else
  {
    l_xmin = xmin;
    l_xmax = xmax;
  }

  /* this function is a refactored version of brent_opt */
  /* if we consider the following structure:
   *
   * 	(a) initialization block
   * 	(b) main loop:
   * 	  (b.1) loop pre-score
   * 	  (b.2) loop score (call to target)
   * 	  (b.3) loop post-score
   *
   * the loop has been removed and it has been split into 2 functions:
   *
   * 1. brent_opt_init, that covers (a) and (b.1)
   * 2. brent_opt_post_loop, that covers (b.3) and (b.1)
   *
   * This way, brent_opt can be refactored as follows:
   *
   * (A) brent_opt_init
   * (B) main loop
   * 	(B.1) loop score (call to target)
   * 	(B.2) brent_opt_post_loop
   *
   * For parallel execution, the synchronization point occurs right at the
   * beginning of the loop.
   *
   * Moreover, the optimization state is encapsulated in a `struct opt_params *`
   * that can be defined as an array holding the optimization state for each
   * local partition.
   *
   * I observed no impact in results nor in runtime
   */

  for (i = 0; i < xnum; ++i)
  {
    double eps;
    int outbounds_ax, outbounds_cx;

    /* first attempt to bracketize minimum */
    if (xguess[i] < l_xmin[i])
      xguess[i] = l_xmin[i];
    if (xguess[i] > l_xmax[i])
      xguess[i] = l_xmax[i];

    // TODO: this is a quick workaround to make it work for xguess==0
    // But we should double-check this bracketing heuristic! (alexey)
    eps = xguess[i] > 0 ? xguess[i] * xtol * 50.0 : 2. * xtol;

    ax[i] = xguess[i] - eps;
    outbounds_ax = ax[i] < l_xmin[i];
    if (outbounds_ax)
      ax[i] = l_xmin[i];
    cx[i] = xguess[i] + eps;
    outbounds_cx = cx[i] > l_xmax[i];
    if (outbounds_cx)
      cx[i] = l_xmax[i];
  }

  target_funk (params, ax, fa, NULL);
  target_funk (params, xguess, fb, NULL);
  target_funk (params, cx, fc, NULL);
  target_funk (params, l_xmin, fxmin, NULL);
  target_funk (params, l_xmax, fxmax, NULL);

  for (i = 0; i < xnum; ++i)
    init_ok[i] = 1;

  /* check if this works */
  for (i = 0; i < xnum; ++i)
  {
    /* skip params that do not need to be optimized */
    if (opt_mask && !opt_mask[i])
      continue;

    /* if it works use these borders else be conservative */
    if ((fa[i] < fb[i]) || (fc[i] < fb[i]))
    {
      fa[i] = fxmin[i];
      fc[i] = fxmax[i];
      ax[i] = l_xmin[i];
      cx[i] = l_xmax[i];
    }

    DBG("param[%d] a x c / fa fx fc: %lf, %lf, %lf / %lf, %lf, %lf / %lf\n",
        i, ax[i], xguess[i], cx[i], fa[i], fb[i], fc[i], l_xmax[i]);

    if (!brent_opt_init(ax[i], xguess[i], cx[i], xtol, fx, f2x,
                        fa[i], fb[i], fc[i], &brent_params[i]))
    {
      init_ok[i] = 0;
      break;
    }
  }

  free(ax);
  free(cx);
  free(fa);
  free(fb);
  free(fc);
  free(fxmin);
  free(fxmax);

  if (global_range)
  {
    free(l_xmin);
    free(l_xmax);
  }

  // check if BRENT initialization has failed in *any* of the threads
  target_funk (params, NULL, NULL, init_ok);
  int init_failed = !init_ok[xnum];
  free(init_ok);

  if (init_failed)
  {
    /* restore the original parameter value */
    target_funk (params, xguess, fx, NULL);
    free(brent_params);
    pllmod_set_error(PLLMOD_OPT_ERROR_BRENT_INIT, "BRENT: initialization failed!");
    return PLL_FAILURE;
  }

  int * converged = (int *) calloc(xnum+1, sizeof(int));
  double * u = (double *) calloc(xnum, sizeof(double));
  double * fu = (double *) calloc(xnum, sizeof(double));

  int iter_num = 0;
  while (iterate)
  {
    for (i = 0; i < xnum; ++i)
      u[i] = brent_params[i].u;

    target_funk (params, u, fu, converged);

//    DBG("iter: %d, u: %lf, fu: %lf\n", iter_num, u[2], fu[2]);

    /* last element in converged[] array is "all converged" flag */
    iterate = !converged[xnum];
    for (i = 0; i < xnum; ++i)
    {
      /* skip params that do not need to be optimized */
      if (opt_mask && !opt_mask[i])
        continue;

      if (!converged[i])
      {
        brent_params[i].fu = fu[i];
        converged[i] = !brent_opt_post_loop(&brent_params[i]);
      }
    }

    iter_num++;
    iterate &= (iter_num <= ITMAX);
  }

  /* if new score is worse, return initial value */
  for (i = 0; i < xnum; ++i)
  {
    xopt[i] = (brent_params[i].fx > brent_params[i].fstartx) ?
        brent_params[i].startx : brent_params[i].x;
  }

  target_funk (params, xopt, fx, NULL);

  DBG("xopt: %lf, LH: %lf\n", xopt[2], fx ? fx[2] : NAN);

  free(u);
  free(fu);
  free(converged);
  free(brent_params);
  return PLL_SUCCESS;
}

// TODO: remove at some point (not used anymore)
#if 0
static double brent_opt (double ax, double bx, double cx, double tol,
                         double *foptx, double *f2optx, double fax,
                         double fbx, double fcx,
                         void * params,
                         double (*target_funk)(
                             void *,
                             double))
{
  int iter;
  double a, b, d = 0, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x,
      xm;
  double xw, wv, vx;
  double e = 0.0;

  a = (ax < cx ? ax : cx);
  b = (ax > cx ? ax : cx);
  x = bx;
  fx = fbx;
  if (fax < fcx)
  {
    w = ax;
    fw = fax;
    v = cx;
    fv = fcx;
  }
  else
  {
    w = cx;
    fw = fcx;
    v = ax;
    fv = fax;
  }

  for (iter = 1; iter <= ITMAX; iter++)
  {
    xm = 0.5 * (a + b);
    tol2 = 2.0 * (tol1 = tol * fabs (x) + ZEPS);
    if (fabs (x - xm) <= (tol2 - 0.5 * (b - a)))
    {
      *foptx = fx;
      xw = x - w;
      wv = w - v;
      vx = v - x;
      *f2optx = 2.0 * (fv * xw + fx * wv + fw * vx)
          / (v * v * xw + x * x * wv + w * w * vx);
      return x;
    }

    if (fabs (e) > tol1)
    {
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = 2.0 * (q - r);
      if (q > 0.0)
        p = -p;
      q = fabs (q);
      etemp = e;
      e = d;
      if (fabs (p) >= fabs (0.5 * q * etemp) || p <= q * (a - x)
          || p >= q * (b - x))
        d = CGOLD * (e = (x >= xm ? a - x : b - x));
      else
      {
        d = p / q;
        u = x + d;
        if (u - a < tol2 || b - u < tol2)
          d = SIGN(tol1, xm - x);
      }
    }
    else
    {
      d = CGOLD * (e = (x >= xm ? a - x : b - x));
    }

    u = (fabs (d) >= tol1 ? x + d : x + SIGN(tol1, d));
    fu = target_funk (params, u);
    if (fu <= fx)
    {
      if (u >= x)
        a = x;
      else
        b = x;

      SHFT(v, w, x, u)
      SHFT(fv, fw, fx, fu)
    }
    else
    {
      if (u < x)
        a = u;
      else
        b = u;
      if (fu <= fw || w == x)
      {
        v = w;
        w = u;
        fv = fw;
        fw = fu;
      }
      else if (fu <= fv || v == x || v == w)
      {
        v = u;
        fv = fu;
      }
    }
  }

  *foptx = fx;
  xw = x - w;
  wv = w - v;
  vx = v - x;
  *f2optx = 2.0 * (fv * xw + fx * wv + fw * vx)
      / (v * v * xw + x * x * wv + w * w * vx);

  return x;
}
#endif

#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT
#undef SIGN

/* most of the code for Brent optimization taken from IQ-Tree
 * http://www.cibiv.at/software/iqtree
 * --------------------------------------------------------------------------
 * Bui Quang Minh, Minh Anh Thi Nguyen, and Arndt von Haeseler (2013)
 * Ultrafast approximation for phylogenetic bootstrap.
 * Mol. Biol. Evol., 30:1188-1195. (free reprint, DOI: 10.1093/molbev/mst024) */

 /**
  * Minimize a single-variable function using Brent.
  *
  * The target function must evaluate the function at a certain point,
  * and it requires 2 parameters: (1) custom data (if needed),
  * and (2) the value where the function is evaluated.
  *
  * @param  xmin       lower bound
  * @param  xguess     first guess for the free variable
  * @param  xmax       upper bound
  * @param  xtol       tolerance of the minimization method
  * @param  params     custom parameters required by the target function
  * @param  deriv_func target function
  *
  * @return            the parameter value that minimizes the function in [xmin,xmax]
  */
PLL_EXPORT double pllmod_opt_minimize_brent(double xmin,
                                           double xguess,
                                           double xmax,
                                           double xtol,
                                           double *fx,
                                           double *f2x,
                                           void * params,
                                           double (*target_funk)(
                                               void *,
                                               double))
{
  double optx = xguess;

  struct brent_wrapper_params wrap_params;
  wrap_params.target_funk = target_funk;
  wrap_params.params = params;

  brent_opt_alt (1, NULL, &xmin, &xguess, &xmax, xtol, &optx, fx, f2x, &wrap_params,
                    target_funk_wrapper, 1);

//  optx = brent_opt(xmin, xguess, xmax, xtol, fx, f2x, params, target_funk);

  return optx; /* return optimal x */
}

/**
 * Run brent optimization for multiple variables in parallel
 * (e.g., unlinked alphas for multiple partitions)
 *
 * @param xnum number of variables/partitions
 * @param opt_mask opt_mask[i]=1/0: optimize/skip parameter i
 * @param xmin minimum value(s) (see @param global_range)
 * @param xguess array of staring values
 * @param xmax maximum value(s) (see @param global_range)
 * @param xtol tolerance
 * @param xopt optimal variable values [out]
 * @param fx auxiliary variable to keep state between runs [out]
 * @param f2x auxiliary variable to keep state between runs [out]
 * @param params parameters to be passed to target_funk
 * @param target_funk target function, parameters: 1) params, 2) array of x values,
 * 3) array of scores [out], 4) convergence flags
 * @param global_range 0=xmin/xmax point to arrays of size xnum with individual
 * per-variable ranges; 1=xmin/xmax is a global range for all variables
 */
PLL_EXPORT int pllmod_opt_minimize_brent_multi(int xnum,
                                               int * opt_mask,
                                               double * xmin,
                                               double * xguess,
                                               double * xmax,
                                               double xtol,
                                               double * xopt,
                                               double * fx,
                                               double * f2x,
                                               void * params,
                                               double (*target_funk)(
                                                       void *,
                                                      double *,
                                                      double *,
                                                      int *),
                                               int global_range)
{
  return brent_opt_alt (xnum, opt_mask, xmin, xguess, xmax, xtol, xopt, fx, f2x, params,
                          target_funk, global_range);
}

/******************************************************************************/
/* EXPECTATION-MAXIMIZATION (EM)     */
/* Wang, Li, Susko, and Roger (2008) */
/******************************************************************************/
PLL_EXPORT void pllmod_opt_minimize_em(double *w,
                                       unsigned int w_count,
                                       double *sitecat_lh,
                                       unsigned int *site_w,
                                       unsigned int l,
                                       void * params,
                                       double (*update_sitecatlk_funk)(
                                           void *,
                                           double *))
{
  unsigned int i, c;
  unsigned int max_steps = 10;
  int converged = 0;
  int ratio_scale = 0;

  double *new_prop = (double *) malloc(sizeof(double) * w_count);
  double *ratio_prop = (double *) malloc(sizeof(double) * w_count);

  while (!converged && max_steps--)
  {
    /* update site-cat LK */
    update_sitecatlk_funk(params, sitecat_lh);

    // Expectation
    double *this_lk_cat = sitecat_lh;
    if (ratio_scale)
    {
      for (i = 0; i < l; ++i)
      {
        for (c = 0; c < w_count; c++)
        {
          this_lk_cat[c] *= ratio_prop[c];
        }
        this_lk_cat += w_count;
      }
    }
    else
      ratio_scale = 1;

    memset (new_prop, 0, w_count * sizeof(double));

    this_lk_cat = sitecat_lh;
    for (i=0; i<l; ++i)
    {
      //TODO: Check for p_invar
      double lk_ptn = 0;
      for (c = 0; c < w_count; c++)
      {
        lk_ptn += this_lk_cat[c];
      }
      lk_ptn = site_w[i] / lk_ptn;
      for (c = 0; c < w_count; c++)
      {
        new_prop[c] += this_lk_cat[c] * lk_ptn;
      }
      this_lk_cat += w_count;
    }

    // Maximization
    converged = 1;
    for (c = 0; c < w_count; c++)
    {
      new_prop[c] /= l;

      // check for convergence
      converged = converged && (fabs (w[c] - new_prop[c]) < 1e-4);
      ratio_prop[c] = new_prop[c] / w[c];
      w[c] = new_prop[c];
    }
  }

  free(ratio_prop);
  free(new_prop);
}
