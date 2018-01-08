/*
 Copyright (C) 2016 Alexey Kozlov, Diego Darriba

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
 * @file models_dna.c
 *
 * @brief
 *
 * @author Alexey Kozlov
 * @author Diego Darriba
 */

#include <string.h>

#include "pllmod_util.h"
#include "../pllmod_common.h"

#define DNA_MODELS_COUNT 22

/*                                 AC AG AT CG CT GT       */
static const double dna_rates_equal[] = {1, 1, 1, 1, 1, 1};

/*                                  A     C     G     T    */
static const double dna_freqs_equal[] = {0.25, 0.25, 0.25, 0.25};

/*                                 A  C  G  T              */
static const int dna_sym_freq_equal[] = {0, 0, 0, 0};
static const int dna_sym_freq_free[]  = {0, 1, 2, 3};

/*                                 AC AG AT CG CT GT       */
static const int dna_sym_rate_equal[] = {0, 0, 0, 0, 0, 0};
static const int dna_sym_rate_free[]  = {0, 1, 2, 3, 4, 5};
static const int dna_sym_rate_tvts[]  = {0, 1, 0, 0, 1, 0};
static const int dna_sym_rate_tn93[]  = {0, 1, 0, 0, 2, 0};
static const int dna_sym_rate_k81[]   = {0, 1, 2, 2, 1, 0};
static const int dna_sym_rate_tpm2[]  = {0, 1, 0, 2, 1, 2};
static const int dna_sym_rate_tpm3[]  = {0, 1, 2, 0, 1, 2};
static const int dna_sym_rate_tim1[]  = {0, 1, 2, 2, 3, 0};
static const int dna_sym_rate_tim2[]  = {0, 1, 0, 2, 3, 2};
static const int dna_sym_rate_tim3[]  = {0, 1, 2, 0, 3, 2};
static const int dna_sym_rate_tvm[]   = {0, 1, 2, 3, 1, 4};

static const pllmod_subst_model_t dna_model_list[DNA_MODELS_COUNT] =
{
    /*       states  model rates         model freqs      rate symmetries     frequencies sym.           */
    {"JC",     4,    dna_rates_equal,   dna_freqs_equal, dna_sym_rate_equal, dna_sym_freq_equal, 0 },

    {"K80",    4,    NULL,              dna_freqs_equal, dna_sym_rate_tvts,  dna_sym_freq_equal, 0 },

    {"F81",    4,    dna_rates_equal,   NULL,            dna_sym_rate_equal, dna_sym_freq_free, 0 },

    {"HKY",    4,    NULL,              NULL,            dna_sym_rate_tvts,  dna_sym_freq_free, 0 },

    {"TN93ef", 4,    NULL,              dna_freqs_equal, dna_sym_rate_tn93,  dna_sym_freq_equal, 0 },

    {"TN93",   4,    NULL,              NULL,            dna_sym_rate_tn93,  dna_sym_freq_free, 0 },

    {"K81",    4,    NULL,              dna_freqs_equal, dna_sym_rate_k81,   dna_sym_freq_equal, 0 },

    {"K81uf",  4,    NULL,              NULL,            dna_sym_rate_k81,   dna_sym_freq_free, 0 },

    {"TPM2",   4,    NULL,              dna_freqs_equal, dna_sym_rate_tpm2,  dna_sym_freq_equal, 0 },

    {"TPM2uf", 4,    NULL,              NULL,            dna_sym_rate_tpm2,  dna_sym_freq_free, 0 },

    {"TPM3",   4,    NULL,              dna_freqs_equal, dna_sym_rate_tpm3,  dna_sym_freq_equal, 0 },

    {"TPM3uf", 4,    NULL,              NULL,            dna_sym_rate_tpm3,  dna_sym_freq_free, 0 },

    {"TIM1",   4,    NULL,              dna_freqs_equal, dna_sym_rate_tim1,  dna_sym_freq_equal, 0 },

    {"TIM1uf", 4,    NULL,              NULL,            dna_sym_rate_tim1,  dna_sym_freq_free, 0 },

    {"TIM2",   4,    NULL,              dna_freqs_equal, dna_sym_rate_tim2,  dna_sym_freq_equal, 0 },

    {"TIM2uf", 4,    NULL,              NULL,            dna_sym_rate_tim2,  dna_sym_freq_free, 0 },

    {"TIM3",   4,    NULL,              dna_freqs_equal, dna_sym_rate_tim3,  dna_sym_freq_equal, 0 },

    {"TIM3uf", 4,    NULL,              NULL,            dna_sym_rate_tim3,  dna_sym_freq_free, 0 },

    {"TVMef",  4,    NULL,              dna_freqs_equal, dna_sym_rate_tvm,   dna_sym_freq_equal, 0 },

    {"TVM",    4,    NULL,              NULL,            dna_sym_rate_tvm,   dna_sym_freq_free, 0 },

    {"SYM",    4,    NULL,              dna_freqs_equal, dna_sym_rate_free,  dna_sym_freq_equal, 0 },

    {"GTR",    4,    NULL,              NULL,            dna_sym_rate_free,  dna_sym_freq_free, 0 }
};

static const pllmod_subst_model_alias_t dna_model_aliases[] =
{
	{"TrNef", "TN93ef"},
	{"TrN", "TN93"},
	{"TPM1", "K81"},
	{"TPM1uf", "K81uf"}
};

const int ALIAS_COUNT =
		sizeof(dna_model_aliases) / sizeof(pllmod_subst_model_alias_t);

static int get_model_index(const char * model_name)
{
  int i;
  const char * resolved_name = model_name;

  /* resolve model aliases first (e.g., TPM1 -> K81) */
  for (i = 0; i < ALIAS_COUNT; ++i)
  {
	if (strcasecmp(model_name, dna_model_aliases[i].alias) == 0)
	{
		resolved_name = dna_model_aliases[i].primary_name;
		break;
	}
  }

  /* search for the model */
  for (i = 0; i < DNA_MODELS_COUNT; ++i)
    if (strcasecmp(resolved_name, dna_model_list[i].name) == 0)
      return i;

  /* model not found*/
  return -1;
}

/**
 * @brief Returns number of available built-in DNA evolution models
 */
PLL_EXPORT unsigned int pllmod_util_model_count_dna()
{
  return DNA_MODELS_COUNT;
}

/**
 * @brief Returns list of available built-in DNA evolution models (names)
 */
PLL_EXPORT char ** pllmod_util_model_names_dna()
{
  char ** names = calloc(DNA_MODELS_COUNT, sizeof(char *));

  int i;
  for (i = 0; i < DNA_MODELS_COUNT; ++i)
    {
      const char * model_name = dna_model_list[i].name;
      names[i] = malloc(strlen(model_name)+1);
      strcpy(names[i], model_name);
    }

  return names;
}

/**
 * @brief Returns 1 if built-in DNA models with a given name exists and 0 otherwise
 */
PLL_EXPORT int pllmod_util_model_exists_dna(const char * model_name)
{
  return get_model_index(model_name) >= 0 ? 1 : 0;
}

/**
 * @brief Returns properties of the specified DNA evolution model
 *
 * See pllmod_model_t definition for details
 *
 * @param model_name name of the DNA model
 *
 * @return model info structure, or NULL if model doesn't exist
 */
PLL_EXPORT pllmod_subst_model_t * pllmod_util_model_info_dna(const char * model_name)
{
  const int model_index = get_model_index(model_name);
  if (model_index >= 0)
    {
      return pllmod_util_model_clone(&dna_model_list[model_index]);
    }
  else
    {
      pllmod_set_error(PLLMOD_UTIL_ERROR_MODEL_UNKNOWN, "DNA model not found: %s", model_name);
      return NULL;
    }
}
