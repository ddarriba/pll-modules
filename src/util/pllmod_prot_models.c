#include <string.h>

#include "pllmod_util.h"
#include "../pllmod_common.h"

#define PROT_MODELS_COUNT 21

static const double * protein_models_rates_list[PROT_MODELS_COUNT] =
 {
   pll_aa_rates_dayhoff,
   pll_aa_rates_lg,
   pll_aa_rates_dcmut,
   pll_aa_rates_jtt,
   pll_aa_rates_mtrev,
   pll_aa_rates_wag,
   pll_aa_rates_rtrev,
   pll_aa_rates_cprev,
   pll_aa_rates_vt,
   pll_aa_rates_blosum62,
   pll_aa_rates_mtmam,
   pll_aa_rates_mtart,
   pll_aa_rates_mtzoa,
   pll_aa_rates_pmb,
   pll_aa_rates_hivb,
   pll_aa_rates_hivw,
   pll_aa_rates_jttdcmut,
   pll_aa_rates_flu,
   pll_aa_rates_stmtrev,
   (double*) pll_aa_rates_lg4m,
   (double*) pll_aa_rates_lg4x
 };

static const double * protein_models_freqs_list[PROT_MODELS_COUNT] =
 {
   pll_aa_freqs_dayhoff,
   pll_aa_freqs_lg,
   pll_aa_freqs_dcmut,
   pll_aa_freqs_jtt,
   pll_aa_freqs_mtrev,
   pll_aa_freqs_wag,
   pll_aa_freqs_rtrev,
   pll_aa_freqs_cprev,
   pll_aa_freqs_vt,
   pll_aa_freqs_blosum62,
   pll_aa_freqs_mtmam,
   pll_aa_freqs_mtart,
   pll_aa_freqs_mtzoa,
   pll_aa_freqs_pmb,
   pll_aa_freqs_hivb,
   pll_aa_freqs_hivw,
   pll_aa_freqs_jttdcmut,
   pll_aa_freqs_flu,
   pll_aa_freqs_stmtrev,
   (double*) pll_aa_freqs_lg4m,
   (double*) pll_aa_freqs_lg4x
 };

static const char * protein_models_names_list[PROT_MODELS_COUNT] =
{
   "DAYHOFF",
   "LG",
   "DCMUT",
   "JTT",
   "MTREV",
   "WAG",
   "RTREV",
   "CPREV",
   "VT",
   "BLOSUM62",
   "MTMAM",
   "MTART",
   "MTZOA",
   "PMB",
   "HIVB",
   "HIVW",
   "JTTDCMUT",
   "FLU",
   "STMTREV",
   "LG4M",
   "LG4X"
};

static const unsigned int protein_models_matrix_count[PROT_MODELS_COUNT] =
{
    /* 19 single-matrix models */
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    4, /* LG4M */
    4  /* LG4X */
};

static int get_model_index(const char * model_name)
{
  int i;
  for (i = 0; i < PROT_MODELS_COUNT; ++i)
    if (strcasecmp(model_name, protein_models_names_list[i]) == 0)
      return i;

  /* model not found*/
  return -1;
}

/**
 * @brief Returns number of available built-in protein evolution models
 */
PLL_EXPORT unsigned int pllmod_util_get_protein_model_count()
{
  return PROT_MODELS_COUNT;
}

/**
 * @brief Returns list of available built-in protein evolution models (names)
 */
PLL_EXPORT const char ** pllmod_util_get_protein_model_names()
{
  return protein_models_names_list;
}

/**
 * @brief Returns 1 if built-in protein models with a given name exists and 0 otherwise
 */
PLL_EXPORT int pllmod_util_protein_model_exists(const char * model_name)
{
  return get_model_index(model_name) >= 0 ? 1 : 0;
}

/**
 * @brief Returns array of model-specific substitution rates
 *
 * Size of this array is 20 * <number of rate matrices> (i.e., 20 for WAG and
 * other single-matrix models, 80 for LG4M and LG4X etc.)
 *
 * @param model_name name of the protein model
 *
 * @return array with AA substitution rates, or NULL if model doesn't exist
 */
PLL_EXPORT const double * pllmod_util_get_protein_model_rates(const char * model_name)
{
  const int model_index = get_model_index(model_name);
  if (model_index >= 0)
    {
      return protein_models_rates_list[model_index];
    }
  else
    {
      pllmod_set_error(PLLMOD_ERROR_UNKNOWN_MODEL, "Protein model not found: %s", model_name);
      return NULL;
    }
}

/**
 * @brief Returns array of model-specific aminoacid frequencies
 *
 * Size of this array is 20 * <number of rate matrices> (i.e., 20 for WAG and
 * other single-matrix models, 80 for LG4M and LG4X etc.)
 *
 * @param model_name name of the protein model
 *
 * @return array of AA frequencies, or NULL if model doesn't exist
 */
PLL_EXPORT const double * pllmod_util_get_protein_model_freqs(const char * model_name)
{
  const int model_index = get_model_index(model_name);
  if (model_index >= 0)
    {
      return protein_models_freqs_list[model_index];
    }
  else
    {
      pllmod_set_error(PLLMOD_ERROR_UNKNOWN_MODEL, "Protein model not found: %s", model_name);
      return NULL;
    }
}

/**
 * @brief Returns number of rate matrices for the model (e.g., LG4X = 4, LG = 1)
 *
 * @param model_name name of the protein model
 *
 * @return number of rate matrices, 0 on error
 */
PLL_EXPORT unsigned int pllmod_util_get_protein_model_matrix_count(const char * model_name)
{
  const int model_index = get_model_index(model_name);
  if (model_index >= 0)
    {
      return protein_models_matrix_count[model_index];
    }
  else
    {
      pllmod_set_error(PLLMOD_ERROR_UNKNOWN_MODEL, "Protein model not found: %s", model_name);
      return 0;
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
PLL_EXPORT int pllmod_util_set_protein_model(pll_partition_t * partition, const char * model_name, int model_freqs)
{
  const int model_index = get_model_index(model_name);
  if (model_index >= 0)
    {
      unsigned int i;
      for (i = 0; i < protein_models_matrix_count[model_index]; ++i)
        {
          pll_set_subst_params(partition, i,
                               protein_models_rates_list[model_index] + i * PLLMOD_PROT_RATES_COUNT);
        }
      if (model_freqs)
        {
          for (i = 0; i < protein_models_matrix_count[model_index]; ++i)
            {
              pll_set_frequencies(partition, i,
                                  protein_models_freqs_list[model_index] + i * PLLMOD_PROT_FREQS_COUNT);
            }
        }
      return PLL_SUCCESS;
    }
  else
    {
      pllmod_set_error(PLLMOD_ERROR_UNKNOWN_MODEL, "Protein model not found: %s", model_name);
      return PLL_FAILURE;
    }
}

