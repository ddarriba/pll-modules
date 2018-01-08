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

/* general single-matrix models */
const pllmod_subst_model_t M_DAYHOFF  = {"DAYHOFF",   20, pll_aa_rates_dayhoff,  pll_aa_freqs_dayhoff,  NULL, NULL, 0 };
const pllmod_subst_model_t M_LG       = {"LG",        20, pll_aa_rates_lg,       pll_aa_freqs_lg,       NULL, NULL, 0 };
const pllmod_subst_model_t M_DCMUT    = {"DCMUT",     20, pll_aa_rates_dcmut,    pll_aa_freqs_dcmut,    NULL, NULL, 0 };
const pllmod_subst_model_t M_JTT      = {"JTT",       20, pll_aa_rates_jtt,      pll_aa_freqs_jtt,      NULL, NULL, 0 };
const pllmod_subst_model_t M_MTREV    = {"MTREV",     20, pll_aa_rates_mtrev,    pll_aa_freqs_mtrev,    NULL, NULL, 0 };
const pllmod_subst_model_t M_WAG      = {"WAG",       20, pll_aa_rates_wag,      pll_aa_freqs_wag,      NULL, NULL, 0 };
const pllmod_subst_model_t M_RTREV    = {"RTREV",     20, pll_aa_rates_rtrev,    pll_aa_freqs_rtrev,    NULL, NULL, 0 };
const pllmod_subst_model_t M_CPREV    = {"CPREV",     20, pll_aa_rates_cprev,    pll_aa_freqs_cprev,    NULL, NULL, 0 };
const pllmod_subst_model_t M_VT       = {"VT",        20, pll_aa_rates_vt,       pll_aa_freqs_vt,       NULL, NULL, 0 };
const pllmod_subst_model_t M_BLOSUM62 = {"BLOSUM62",  20, pll_aa_rates_blosum62, pll_aa_freqs_blosum62, NULL, NULL, 0 };
const pllmod_subst_model_t M_MTMAM    = {"MTMAM",     20, pll_aa_rates_mtmam,    pll_aa_freqs_mtmam,    NULL, NULL, 0 };
const pllmod_subst_model_t M_MTART    = {"MTART",     20, pll_aa_rates_mtart,    pll_aa_freqs_mtart,    NULL, NULL, 0 };
const pllmod_subst_model_t M_MTZOA    = {"MTZOA",     20, pll_aa_rates_mtzoa,    pll_aa_freqs_mtzoa,    NULL, NULL, 0 };
const pllmod_subst_model_t M_PMB      = {"PMB",       20, pll_aa_rates_pmb,      pll_aa_freqs_pmb,      NULL, NULL, 0 };
const pllmod_subst_model_t M_HIVB     = {"HIVB",      20, pll_aa_rates_hivb,     pll_aa_freqs_hivb,     NULL, NULL, 0 };
const pllmod_subst_model_t M_HIVW     = {"HIVW",      20, pll_aa_rates_hivw,     pll_aa_freqs_hivw,     NULL, NULL, 0 };
const pllmod_subst_model_t M_JTTDCMUT = {"JTT-DCMUT", 20, pll_aa_rates_jttdcmut, pll_aa_freqs_jttdcmut, NULL, NULL, 0 };
const pllmod_subst_model_t M_FLU      = {"FLU",       20, pll_aa_rates_flu,      pll_aa_freqs_flu,      NULL, NULL, 0 };
const pllmod_subst_model_t M_STMTREV  = {"STMTREV",   20, pll_aa_rates_stmtrev,  pll_aa_freqs_stmtrev,  NULL, NULL, 0 };

/* LG4M components */
const pllmod_subst_model_t M_LG4M1    = {"LG4M1",     20, pll_aa_rates_lg4m[0],  pll_aa_freqs_lg4m[0],  NULL, NULL, 0 };
const pllmod_subst_model_t M_LG4M2    = {"LG4M2",     20, pll_aa_rates_lg4m[1],  pll_aa_freqs_lg4m[1],  NULL, NULL, 0 };
const pllmod_subst_model_t M_LG4M3    = {"LG4M3",     20, pll_aa_rates_lg4m[2],  pll_aa_freqs_lg4m[2],  NULL, NULL, 0 };
const pllmod_subst_model_t M_LG4M4    = {"LG4M4",     20, pll_aa_rates_lg4m[3],  pll_aa_freqs_lg4m[3],  NULL, NULL, 0 };

/* LG4X components */
const pllmod_subst_model_t M_LG4X1    = {"LG4X1",     20, pll_aa_rates_lg4x[0],  pll_aa_freqs_lg4x[0],  NULL, NULL, 0 };
const pllmod_subst_model_t M_LG4X2    = {"LG4X2",     20, pll_aa_rates_lg4x[1],  pll_aa_freqs_lg4x[1],  NULL, NULL, 0 };
const pllmod_subst_model_t M_LG4X3    = {"LG4X3",     20, pll_aa_rates_lg4x[2],  pll_aa_freqs_lg4x[2],  NULL, NULL, 0 };
const pllmod_subst_model_t M_LG4X4    = {"LG4X4",     20, pll_aa_rates_lg4x[3],  pll_aa_freqs_lg4x[3],  NULL, NULL, 0 };

const pllmod_subst_model_t M_PROTGTR = {"PROTGTR",    20, NULL, NULL, NULL, NULL, 0 };


static const pllmod_subst_model_t * prot_model_list[] =
  {
    &M_DAYHOFF,
    &M_LG,
    &M_DCMUT,
    &M_JTT,
    &M_MTREV,
    &M_WAG,
    &M_RTREV,
    &M_CPREV,
    &M_VT,
    &M_BLOSUM62,
    &M_MTMAM,
    &M_MTART,
    &M_MTZOA,
    &M_PMB,
    &M_HIVB,
    &M_HIVW,
    &M_JTTDCMUT,
    &M_FLU,
    &M_STMTREV,

    &M_LG4M1,
    &M_LG4M2,
    &M_LG4M3,
    &M_LG4M4,

    &M_LG4X1,
    &M_LG4X2,
    &M_LG4X3,
    &M_LG4X4,

    &M_PROTGTR
  };

const int PROT_MODELS_COUNT = sizeof(prot_model_list) / sizeof(pllmod_subst_model_t *);

/* mixture models */
const pllmod_subst_model_t * lg4m_matrices[] = {&M_LG4M1, &M_LG4M2, &M_LG4M3, &M_LG4M4};
const pllmod_mixture_model_t M_LG4M = {"LG4M", 4, (pllmod_subst_model_t **) lg4m_matrices,
    NULL, NULL, PLLMOD_UTIL_MIXTYPE_GAMMA };

const pllmod_subst_model_t * lg4x_matrices[] = {&M_LG4X1, &M_LG4X2, &M_LG4X3, &M_LG4X4};
const pllmod_mixture_model_t M_LG4X = {"LG4X", 4, (pllmod_subst_model_t **) lg4x_matrices,
    NULL, NULL, PLLMOD_UTIL_MIXTYPE_FREE };

static const pllmod_mixture_model_t * protmix_model_list[] =
  {
      &M_LG4M,
      &M_LG4X
  };

const int PROTMIX_MODELS_COUNT = sizeof(protmix_model_list) / sizeof(pllmod_mixture_model_t *);


static int get_model_index(const char * model_name)
{
  int i;
  for (i = 0; i < PROT_MODELS_COUNT; ++i)
    if (strcasecmp(model_name, prot_model_list[i]->name) == 0)
      return i;

  /* model not found*/
  return -1;
}

static int get_mixmodel_index(const char * model_name)
{
  int i;
  for (i = 0; i < PROTMIX_MODELS_COUNT; ++i)
    if (strcasecmp(model_name, protmix_model_list[i]->name) == 0)
      return i;

  /* model not found*/
  return -1;
}

/**
 * @brief Returns number of available built-in protein evolution models
 */
PLL_EXPORT unsigned int pllmod_util_model_count_protein()
{
  return PROT_MODELS_COUNT;
}

/**
 * @brief Returns list of available built-in protein evolution models (names)
 */
PLL_EXPORT char ** pllmod_util_model_names_protein()
{
  char ** names = calloc(PROT_MODELS_COUNT, sizeof(char *));

  int i;
  for (i = 0; i < PROT_MODELS_COUNT; ++i)
    {
      const char * model_name = prot_model_list[i]->name;
      names[i] = malloc(strlen(model_name)+1);
      strcpy(names[i], model_name);
    }

  return names;
}

/**
 * @brief Returns 1 if built-in protein models with a given name exists and 0 otherwise
 */
PLL_EXPORT int pllmod_util_model_exists_protein(const char * model_name)
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
PLL_EXPORT pllmod_subst_model_t * pllmod_util_model_info_protein(const char * model_name)
{
  const int model_index = get_model_index(model_name);
  if (model_index >= 0)
    {
      return pllmod_util_model_clone(prot_model_list[model_index]);
    }
  else
    {
      pllmod_set_error(PLLMOD_UTIL_ERROR_MODEL_UNKNOWN, "Protein model not found: %s", model_name);
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
PLL_EXPORT int pllmod_util_model_set_protein(pll_partition_t * partition, const char * model_name, int model_freqs)
{
  const int model_index = get_model_index(model_name);
  if (model_index >= 0)
    {
      pll_set_subst_params(partition, 0, prot_model_list[model_index]->rates);
      if (model_freqs)
        {
          pll_set_frequencies(partition, 0, prot_model_list[model_index]->freqs);
        }
      return PLL_SUCCESS;
    }
  else
    {
      pllmod_set_error(PLLMOD_UTIL_ERROR_MODEL_UNKNOWN, "Protein model not found: %s", model_name);
      return PLL_FAILURE;
    }
}

PLL_EXPORT int pllmod_util_model_exists_protmix(const char * model_name)
{
  return get_mixmodel_index(model_name) >= 0 ? 1 : 0;
}

PLL_EXPORT pllmod_mixture_model_t * pllmod_util_model_info_protmix(const char * model_name)
{
  const int model_index = get_mixmodel_index(model_name);
  if (model_index >= 0)
    {
      return pllmod_util_model_mixture_clone(protmix_model_list[model_index]);
    }
  else
    {
      pllmod_set_error(PLLMOD_UTIL_ERROR_MODEL_UNKNOWN, "Protein mixture model not found: %s", model_name);
      return NULL;
    }
}

/**
 * @brief Set given protein mixture model to a pll_partition_t instance
 *
 * @param partition partition instance
 * @param model_name name of the protein mixture model
 * @param model_freqs 0: set model rate matrices only, 1: set model AA frequencies as well
 *
 * @return PLL_SUCCESS on success, PLL_FAILURE on error (check pll_errmsg for details)
 */
PLL_EXPORT int pllmod_util_model_set_protmix(pll_partition_t * partition,
                                             const char * model_name, int model_freqs)
{
  const int model_index = get_model_index(model_name);
  if (model_index >= 0)
    {
      pllmod_mixture_model_t * mixture = NULL;

      if (partition->rate_matrices != mixture->ncomp)
        {
          pllmod_set_error(PLLMOD_UTIL_ERROR_MIXTURE_INVALID_SIZE,
                           "Number of partition matrices (%d) differs "
                           "from the number of mixture components (%d)",
                           partition->rate_matrices, mixture->ncomp);
          return PLL_FAILURE;
        }

      unsigned int i;
      for (i = 0; i < mixture->ncomp; ++i)
        {
          pll_set_subst_params(partition, i, mixture->models[i]->rates);
        }
      if (model_freqs)
        {
          for (i = 0; i < mixture->ncomp; ++i)
            {
              pll_set_frequencies(partition, i, mixture->models[i]->freqs);
            }
        }
      return PLL_SUCCESS;
    }
  else
    {
      pllmod_set_error(PLLMOD_UTIL_ERROR_MODEL_UNKNOWN,
                       "Protein model not found: %s", model_name);
      return PLL_FAILURE;
    }
}
