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
#include <string.h>

#include "pllmod_util.h"
#include "../pllmod_common.h"

/*                                      AA CC GG TT AC AG CG AT CT GT          */
static const double gt_rates_equal[] = {   0, 0, 0, 1, 1, 0, 1, 0, 0,    /* AA */
                                              0, 0, 1, 0, 1, 0, 1, 0,    /* CC */
                                                 0, 0, 1, 1, 0, 0, 1,    /* GG */
                                                    0, 0, 0, 1, 1, 1,    /* TT */
                                                       1, 1, 1, 1, 0,    /* AC */
                                                          1, 1, 0, 1,    /* AG */
                                                             0, 1, 1,    /* CG */
                                                                1, 1,    /* AT */
                                                                   1 };  /* CT */


/*                                       A    C    G    T    4    5    6    7    8    9  */
static const double gt_freqs_equal[] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 };

/*                                 A  C  G  T              */
//static int gt_sym_freq_equal[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
//static int gt_sym_freq_free[]  = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

/*                               AA CC GG TT AC AG CG AT CT GT          */
static int gt_sym_rate_free1[] = {   0, 0, 0, 1, 2, 0, 3, 0, 0,    /* AA */
                                        0, 0, 4, 0, 5, 0, 6, 0,    /* CC */
                                           0, 0, 7, 8, 0, 0, 9,    /* GG */
                                              0, 0, 0, 9,10,11,    /* TT */
                                                12,13,14,15, 0,    /* AC */
                                                   16,17, 0,18,    /* AG */
                                                       0,19,20,    /* CG */
                                                         21,22,    /* AT */
                                                            23 };  /* CT */

/*                                           states  model rates         model freqs   rate symmetries   freq. sym.           */
const pllmod_subst_model_t M_GTM0   = {"GTM0",   10, gt_rates_equal,  NULL,            NULL,              NULL };
const pllmod_subst_model_t M_GTGTR1 = {"GTGTR1", 10, NULL,            NULL,            gt_sym_rate_free1, NULL };
const pllmod_subst_model_t M_GTGTR  = {"GTGTR",  10, NULL,            NULL,            NULL,              NULL };

static const pllmod_subst_model_t * gt_model_list[] =
  {
    &M_GTM0,
    &M_GTGTR1,
    &M_GTGTR
  };

const int GT_MODELS_COUNT = sizeof(gt_model_list) / sizeof(pllmod_subst_model_t *);


static int get_model_index(const char * model_name)
{
  int i;
  for (i = 0; i < GT_MODELS_COUNT; ++i)
    if (strcasecmp(model_name, gt_model_list[i]->name) == 0)
      return i;

  /* model not found*/
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
      const char * model_name = gt_model_list[i]->name;
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
      return pllmod_util_model_clone(gt_model_list[model_index]);
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
      pll_set_subst_params(partition, 0, gt_model_list[model_index]->rates);
      if (model_freqs)
        {
          pll_set_frequencies(partition, 0, gt_model_list[model_index]->freqs);
        }
      return PLL_SUCCESS;
    }
  else
    {
      pllmod_set_error(PLLMOD_UTIL_ERROR_MODEL_UNKNOWN, "Genotype model not found: %s", model_name);
      return PLL_FAILURE;
    }
}
