/*
 Copyright (C) 2016 Alexey Kozlov

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
 * @file models_gt.c
 *
 * @brief
 *
 * @author Alexey Kozlov
 */

#include <string.h>

#include "pllmod_util.h"
#include "../pllmod_common.h"

/*                                       AA CC GG TT AC AG AT CG CT GT          */
static const double gt_rates_equal_sm[] = {  0, 0, 0, 1, 1, 1, 0, 0, 0,    /* AA */
                                                0, 0, 1, 0, 0, 1, 1, 0,    /* CC */
                                                   0, 0, 1, 0, 1, 0, 1,    /* GG */
                                                      0, 0, 1, 0, 1, 1,    /* TT */
                                                         1, 1, 1, 1, 0,    /* AC */
                                                            1, 1, 0, 1,    /* AG */
                                                               0, 1, 1,    /* AT */
                                                                  1, 1,    /* CG */
                                                                     1 };  /* CT */

/*                                      AA CC GG TT AC AG AT CG CT GT          */
static const double gt_rates_equal[] =  {   1, 1, 1, 1, 1, 1, 1, 1, 1,    /* AA */
                                               1, 1, 1, 1, 1, 1, 1, 1,    /* CC */
                                                  1, 1, 1, 1, 1, 1, 1,    /* GG */
                                                     1, 1, 1, 1, 1, 1,    /* TT */
                                                        1, 1, 1, 1, 1,    /* AC */
                                                           1, 1, 1, 1,    /* AG */
                                                              1, 1, 1,    /* AT */
                                                                 1, 1,    /* CG */
                                                                    1 };  /* CT */

/*                                    AA CC GG TT AC AG AT CG CT GT CA GA TA GC TC TG          */
static const double gt16_rates_equal[120] =
                                    {     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,   /* AA */
                                             1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,   /* CC */
                                                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,   /* GG */
                                                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,   /* TT */
                                                      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,   /* AC */
                                                         1, 1, 1, 1, 1, 1, 1, 1, 1, 1,   /* AG */
                                                            1, 1, 1, 1, 1, 1, 1, 1, 1,   /* AT */
                                                               1, 1, 1, 1, 1, 1, 1, 1,   /* CG */
                                                                  1, 1, 1, 1, 1, 1, 1,   /* CT */
                                                                     1, 1, 1, 1, 1, 1,   /* GT */
                                                                        1, 1, 1, 1, 1,   /* CA */
                                                                           1, 1, 1, 1,   /* GA */
                                                                              1, 1, 1,   /* TA */
                                                                                 1, 1,   /* GC */
                                                                                    1};  /* TC */


/*                                      AA   CC   GG   TT   AC   AG   AT   CG   CT   GT */
static const double gt_freqs_equal[] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 };

#define ONE_16 1./16
static const double gt16_freqs_equal[16] = {ONE_16, ONE_16, ONE_16, ONE_16, ONE_16, ONE_16,
                                            ONE_16, ONE_16, ONE_16, ONE_16, ONE_16, ONE_16,
                                            ONE_16, ONE_16, ONE_16, ONE_16};

/*                                 A  C  G  T              */
//static int gt_sym_freq_equal[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
//static int gt_sym_freq_free[]  = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

/*                                  AA  CC  GG  TT  AC  AG  AT  CG  CT  GT         */
static int gt_sym_rate_free_sm[] = {    -1, -1, -1,  0,  1,  2, -1, -1, -1,    /* AA */
                                            -1, -1,  3, -1, -1,  4,  5, -1,    /* CC */
                                                -1, -1,  6, -1,  7, -1,  8,    /* GG */
                                                    -1, -1,  9, -1, 10, 11,    /* TT */
                                                        12, 13, 14, 15, -1,    /* AC */
                                                            16, 17, -1, 18,    /* AG */
                                                                -1, 19, 20,    /* AT */
                                                                    21, 22,    /* CG */
                                                                        23 };  /* CT */

/* A-C: 0, A-G: 1, A-T: 2, C-G: 3, C-T: 4, G-T: 5, others: -1 (not optimize) */
/*                               AA  CC  GG  TT  AC  AG  AT  CG  CT  GT         */
static int gt_sym_rate_dna4[] = {    -1, -1, -1,  1,  1,  2, -1, -1, -1,  /* AA */
                                         -1, -1,  1, -1, -1,  3,  4, -1,  /* CC */
                                             -1, -1,  1, -1,  3, -1,  5,  /* GG */
                                                 -1, -1,  2, -1,  4,  5,  /* TT */
                                                      3,  4,  1,  2, -1,  /* AC */
                                                          5,  0, -1,  2,  /* AG */
                                                             -1,  0,  1,  /* AT */
                                                                  5,  4,  /* CG */
                                                                      3}; /* CT */


/* A-C = A-T = C-G = G-T: 0, A-G = C-T: 1, others: -1 (not optimize) */
/*                               AA  CC  GG  TT  AC  AG  AT  CG  CT  GT         */
static int gt_sym_rate_hky4[] =  {   -1, -1, -1, 0,  1,  0, -1, -1, -1,    /* AA */
                                         -1, -1, 0, -1, -1,  0,  1, -1,    /* CC */
                                            -1, -1,  1, -1,  0, -1,  0,    /* GG */
                                                -1, -1,  0, -1,  1,  1,    /* TT */
                                                     0,  1,  1,  0, -1,    /* AC */
                                                         0,  0, -1,  0,    /* AG */
                                                            -1,  0,  1,    /* AT */
                                                                 0,  1,    /* CG */
                                                                     0 };  /* CT */


/* A-C: 0, A-G: 1, A-T: 2, C-G: 3, C-T: 4, G-T: 5, others: -1 (not optimize) */
/*                                    AA  CC  GG  TT  AC  AG  AT  CG  CT  GT  CA  GA  TA  GC  TC  TG          */
static int gt16_sym_rate_dna4[] =  {      -1, -1, -1,  0,  1,  2, -1, -1, -1,  0,  1,  2, -1, -1, -1,   /* AA */
                                              -1, -1,  0, -1, -1,  3,  4, -1,  0, -1, -1,  3,  4, -1,   /* CC */
                                                  -1, -1,  1, -1,  3, -1,  5, -1,  1, -1,  3, -1,  5,   /* GG */
                                                      -1, -1,  2, -1,  4,  5, -1, -1,  2, -1,  4,  5,   /* TT */
                                                           3,  4,  1,  2, -1, -1, -1, -1,  1,  2, -1,   /* AC */
                                                               5,  0, -1,  2, -1, -1, -1, -1, -1,  2,   /* AG */
                                                                  -1,  0,  1, -1, -1, -1, -1, -1, -1,   /* AT */
                                                                       5,  4,  1, -1, -1, -1, -1,  4,   /* CG */
                                                                           3,  2, -1, -1, -1, -1, -1,   /* CT */
                                                                              -1,  2, -1,  4, -1, -1,   /* GT */
                                                                                   3,  4, -1, -1, -1,   /* CA */
                                                                                       5,  0, -1, -1,   /* GA */
                                                                                          -1,  0,  1,   /* TA */
                                                                                               5, -1,   /* GC */
                                                                                                   3    /* TC */
};


static const pllmod_subst_model_t gt_model_list[] =
{
/*  name    states  model rates         model freqs   rate symmetries   freq. sym.           */
  {"GT10",       10, NULL,               NULL,              gt_sym_rate_dna4,    NULL, 0 },
  {"GT10JC-SM",  10, gt_rates_equal_sm,  gt_freqs_equal,    NULL,                NULL, 0 },
  {"GT10JC",     10, gt_rates_equal,     gt_freqs_equal,    NULL,                NULL, 0 },
  {"GT10GTR-SM", 10, NULL,               NULL,              gt_sym_rate_free_sm, NULL, 0 },
  {"GT10HKY",    10, NULL,               NULL,              gt_sym_rate_hky4,    NULL, 0 },
  {"GT10GTR",    10, NULL,               NULL,              NULL,                NULL, 0 },
  {"GT16",       16, NULL,               NULL,              gt16_sym_rate_dna4,  NULL, 0 },
  {"GT16JC",     16, gt16_rates_equal,   gt16_freqs_equal,  NULL,                NULL, 0 },
  {"GT16GTR",    16, NULL,               NULL,              NULL,                NULL, 0 }
};

const int GT_MODELS_COUNT = sizeof(gt_model_list) / sizeof(pllmod_subst_model_t);

static const pllmod_subst_model_alias_t gt_model_aliases[] =
{
  {"GTJC",     "GT10JC"},
  {"GTJC-SM",  "GT10JC-SM"},
  {"GTGTR4",   "GT10"},
  {"GTGTR",    "GT10GTR"},
  {"GTGTR-SM", "GT10GTR-SM"},
  {"GTHKY4",   "GT10HKY"},
  {"GPGTR4",   "GT16"}
};

static const int ALIAS_COUNT =
                sizeof(gt_model_aliases) / sizeof(pllmod_subst_model_alias_t);


static int get_model_index(const char * model_name)
{
  int i;
  const char * resolved_name = model_name;

  /* resolve model aliases first (e.g., GTGTR4 -> GT10) */
  for (i = 0; i < ALIAS_COUNT; ++i)
  {
    if (strcasecmp(model_name, gt_model_aliases[i].alias) == 0)
    {
      resolved_name = gt_model_aliases[i].primary_name;
      break;
    }
  }

  /* search for the model */
  for (i = 0; i < GT_MODELS_COUNT; ++i)
    if (strcasecmp(resolved_name, gt_model_list[i].name) == 0)
      return i;

  /* model not found */
  return -1;
}

/**
 * @brief Returns number of available built-in protein evolution models
 */
PLL_EXPORT unsigned int pllmod_util_model_count_genotype()
{
  return GT_MODELS_COUNT;
}

/**
 * @brief Returns list of available built-in protein evolution models (names)
 */
PLL_EXPORT char ** pllmod_util_model_names_genotype()
{
  char ** names = calloc(GT_MODELS_COUNT, sizeof(char *));

  int i;
  for (i = 0; i < GT_MODELS_COUNT; ++i)
    {
      const char * model_name = gt_model_list[i].name;
      names[i] = malloc(strlen(model_name)+1);
      strcpy(names[i], model_name);
    }

  return names;
}

/**
 * @brief Returns 1 if built-in genotype models with a given name exists and 0 otherwise
 */
PLL_EXPORT int pllmod_util_model_exists_genotype(const char * model_name)
{
  return get_model_index(model_name) >= 0 ? 1 : 0;
}

PLL_EXPORT int pllmod_util_model_exists_genotype10(const char * model_name)
{
  int idx = get_model_index(model_name);
  if (idx >= 0)
    return gt_model_list[idx].states == 10 ? 1 : 0;
  else
    return 0;
}

PLL_EXPORT int pllmod_util_model_exists_genotype16(const char * model_name)
{
  int idx = get_model_index(model_name);
  if (idx >= 0)
    return gt_model_list[idx].states == 16 ? 1 : 0;
  else
    return 0;
}

/**
 * @brief Returns properties of the specified AA evolution model
 *
 * See pllmod_subst_model_t definition for details
 *
 * @param model_name name of the AA model
 *
 * @return model info structure, or NULL if model doesn't exist
 */
PLL_EXPORT pllmod_subst_model_t * pllmod_util_model_info_genotype(const char * model_name)
{
  const int model_index = get_model_index(model_name);
  if (model_index >= 0)
    {
      return pllmod_util_model_clone(&gt_model_list[model_index]);
    }
  else
    {
      pllmod_set_error(PLLMOD_UTIL_ERROR_MODEL_UNKNOWN, "Genotype model not found: %s", model_name);
      return NULL;
    }
}

/**
 * @brief Set given protein model to a pll_partition_t instance
 *
 * @param partition partition instance
 * @param model_name name of the protein model
 * @param model_freqs 0: set model rate matrices only, 1: set model AA frequencies as well
 *
 * @return PLL_SUCCESS on success, PLL_FAILURE on error (check pll_errmsg for details)
 */
PLL_EXPORT int pllmod_util_model_set_genotype(pll_partition_t * partition, const char * model_name, int model_freqs)
{
  const int model_index = get_model_index(model_name);
  if (model_index >= 0)
    {
      pll_set_subst_params(partition, 0, gt_model_list[model_index].rates);
      if (model_freqs)
        {
          pll_set_frequencies(partition, 0, gt_model_list[model_index].freqs);
        }
      return PLL_SUCCESS;
    }
  else
    {
      pllmod_set_error(PLLMOD_UTIL_ERROR_MODEL_UNKNOWN, "Genotype model not found: %s", model_name);
      return PLL_FAILURE;
    }
}
