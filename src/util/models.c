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

void const_free(const void* ptr)
{
  free((void*) ptr);
}

int * clone_int_array(const int * src, size_t len)
{
  int * dst = (int *) malloc(len * sizeof(int));
  memcpy(dst, src, len * sizeof(int));
  return dst;
}

double * clone_double_array(const double * src, size_t len)
{
  double * dst = (double *) malloc(len * sizeof(double));
  memcpy(dst, src, len * sizeof(double));
  return dst;
}

/* Returns the number of substitution rates for a given number of states */
PLL_EXPORT unsigned int pllmod_util_subst_rate_count(unsigned int states)
{
  return states * (states - 1) / 2;
}

double * pllmod_util_get_equal_freqs(unsigned int states)
{
  double * basefreqs = calloc(states, sizeof(double));
  if (!basefreqs)
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC, "Cannot allocate memory.");
    return NULL;
  }

  unsigned int i;
  for (i = 0; i < states; ++i)
    basefreqs[i] = 1. / states;
  return basefreqs;
}

double * pllmod_util_get_equal_rates(unsigned int states)
{
  const unsigned int rates = pllmod_util_subst_rate_count(states);
  double * substrates = calloc(rates, sizeof(double));
  if (!substrates)
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC, "Cannot allocate memory.");
    return NULL;
  }

  unsigned int i;
  for (i = 0; i < rates; ++i)
    substrates[i] = 1.;
  return substrates;
}


/* @brief Converts string representation of rate symmetries into an array of indices
 *
 *  Identical chars in the input string encode linked rates.
 *  Indices in the output array are 0-based.
 *
 *  Examples for DNA:
 *
 *  "121121" -> [0, 1, 0, 0, 1, 0]  (K80)
 *  "abcdef" -> [0, 1, 2, 3, 4, 5]  (GTR)
 *  "XXXYYY" -> [0, 0, 0, 1, 1, 1]  (user model)
 *
 * @return array defining substitution rate symmetries
 */
PLL_EXPORT int * pllmod_util_model_string_to_sym(const char * s)
{
  size_t len = strlen(s);
  int * sym_list = calloc(len, sizeof(int));
  size_t i;

  int min = s[0];
  for (i = 1; i < len; ++i)
    if (s[i] < min)
      min = s[i];

  for (i = 0; i < len; ++i)
    sym_list[i] = s[i] - min;

  return sym_list;
}

/**
 * @brief Creates a custom substitution model with given parameters
 *
 * @param name name of the custom model
 * @param states number of states (e.g., 4 for DNA)
 * @param rates model substitution rates (NULL=undefined/estimated)
 * @param freqs model base frequencies (NULL=undefined/estimated)
 * @param rate_sym_str string coding substitution rate symmetries
 * @param freq_sym_str string coding base frequencies symmetries
 *
 * @return custom model instance
 */
PLL_EXPORT pllmod_subst_model_t * pllmod_util_model_create_custom(const char * name,
                                                                  unsigned int states,
                                                                  const double * rates,
                                                                  const double * freqs,
                                                                  const char * rate_sym_str,
                                                                  const char * freq_sym_str)
{
  if (states <= 1)
  {
    pllmod_set_error(PLLMOD_UTIL_ERROR_MODEL_INVALID_DEF,
                     "Invalid number of states: %d", states);
    return NULL;
  }

  const size_t rate_count = states * (states - 1) / 2;
  if (rate_sym_str && strlen(rate_sym_str) != rate_count)
  {
    pllmod_set_error(PLLMOD_UTIL_ERROR_MODEL_INVALID_DEF,
                     "Invalid rates symmetry definition: %s", rate_sym_str);
    return NULL;
  }

  if (freq_sym_str && strlen(freq_sym_str) != states)
  {
    pllmod_set_error(PLLMOD_UTIL_ERROR_MODEL_INVALID_DEF,
                     "Invalid freqs symmetry definition: %s", freq_sym_str);
    return NULL;
  }

  pllmod_subst_model_t * model = calloc(1, sizeof(pllmod_subst_model_t));

  model->states = states;
  model->dynamic_malloc = 1;

  if (name)
    model->name = strdup(name);

  model->rates = rates;
  model->freqs = freqs;

  if (rate_sym_str)
    model->rate_sym = pllmod_util_model_string_to_sym(rate_sym_str);

  if (freq_sym_str)
    model->freq_sym = pllmod_util_model_string_to_sym(freq_sym_str);

  return model;
}

/**
 * @brief Creates a copy of substitution model instance
*/
PLL_EXPORT pllmod_subst_model_t * pllmod_util_model_clone(const pllmod_subst_model_t * src)
{
  if (!src)
    return NULL;

  const size_t rate_count = src->states * (src->states - 1) / 2;

  pllmod_subst_model_t * dst = calloc(1, sizeof(pllmod_subst_model_t));

  dst->dynamic_malloc = 1;
  dst->states = src->states;

  if (src->name)
    dst->name = strdup(src->name);

  if (src->rates)
    dst->rates = clone_double_array(src->rates, rate_count);

  if (src->freqs)
    dst->freqs = clone_double_array(src->freqs, src->states);

  if (src->rate_sym)
    dst->rate_sym = clone_int_array(src->rate_sym, rate_count);

  if (src->freq_sym)
    dst->freq_sym = clone_int_array(src->freq_sym, src->states);

  return dst;
}

/**
 * @brief Destroy a substitution model instance and free associated memory
*/
PLL_EXPORT void pllmod_util_model_destroy(pllmod_subst_model_t * model)
{
  if (model->dynamic_malloc)
  {
    if (model->name)
      const_free(model->name);

    if (model->rates)
      const_free(model->rates);

    if (model->freqs)
      const_free(model->freqs);

    if (model->rate_sym)
      const_free(model->rate_sym);

    if (model->freq_sym)
      const_free(model->freq_sym);

    const_free(model);
  }
}

/**
 * @brief Creates a mixture model with given parameters
 *
 * @param name name of the mixture model
 * @param ncomp number of mixture components or matrices (e.g., 4 for LG4X)
 * @param models array of size ncomp with pointers to the mixture components
 * @param mix_rates per-component rates (NULL=undefined/estimated)
 * @param mix_weights per-component weights (NULL=undefined/estimated)
 * @param mix_type defines how rates & weights are estimated:
 *        PLLMOD_MIXTYPE_FIXED = fixed user-specified rates and weights
 *        PLLMOD_MIXTYPE_GAMMA = GAMMA distribution of rates. equal weights
 *        PLLMOD_MIXTYPE_FREE  = rates and weights are estimated by ML
 *
 * @return mixture model instance
 */
PLL_EXPORT pllmod_mixture_model_t * pllmod_util_model_mixture_create(const char * name,
                                                                     unsigned int ncomp,
                                                                     pllmod_subst_model_t ** const models,
                                                                     const double * mix_rates,
                                                                     const double * mix_weights,
                                                                     int mix_type)
{
  if (ncomp <= 0)
  {
    pllmod_set_error(PLLMOD_UTIL_ERROR_MIXTURE_INVALID_SIZE, "Invalid number of components: %d", ncomp);
    return NULL;
  }

  /* check that all components have the same number of states */
  size_t i;
  for (i = 0; i < ncomp; ++i)
  {
    if (models[i]->states != models[0]->states)
      {
        pllmod_set_error(PLLMOD_UTIL_ERROR_MIXTURE_INVALID_COMPONENT,
                         "Distinct number of states in mixture components 0 and %d: %d != %d",
                         i, models[0]->states, models[i]->states);
        return NULL;
      }
  }

  pllmod_mixture_model_t * mixture = calloc(1, sizeof(pllmod_mixture_model_t));

  mixture->ncomp = ncomp;
  mixture->mix_type = mix_type;

  mixture->models = calloc(ncomp, sizeof(pllmod_subst_model_t *));

  for (i = 0; i < ncomp; ++i)
    mixture->models[i] = models[i]->dynamic_malloc ? pllmod_util_model_clone(models[i]) : models[i];

  if (name)
    mixture->name = strdup(name);

  if (mix_rates)
    mixture->mix_rates = clone_double_array(mix_rates, ncomp);

  if (mix_weights)
    mixture->mix_weights = clone_double_array(mix_weights, ncomp);

  return mixture;
}

/**
 * @brief Create a copy of mixture model
*/
PLL_EXPORT pllmod_mixture_model_t * pllmod_util_model_mixture_clone(const pllmod_mixture_model_t * src)
{
  if (src)
    return pllmod_util_model_mixture_create(src->name, src->ncomp, src->models,
                                            src->mix_rates, src->mix_weights, src->mix_type);
  else
    return NULL;
}

/**
 * @brief Destroy a mixture model instance and free associated memory
*/
PLL_EXPORT void pllmod_util_model_mixture_destroy(pllmod_mixture_model_t * mixture)
{
  if (mixture->name)
    free(mixture->name);

  if (mixture->mix_rates)
    free(mixture->mix_rates);

  if (mixture->mix_weights)
    free(mixture->mix_weights);

  if (mixture->models)
  {
    size_t i;
    for (i = 0; i < mixture->ncomp; ++i)
      pllmod_util_model_destroy(mixture->models[i]);

    free(mixture->models);
  }

  free(mixture);
}

/**
 * @brief Creates a custom libpll character map (ASCII code -> bit-encoded state)
 *
 * @param states number of states
 * @param statechars characters that encode states, in respective order (e.g. "ACGT")
 * @param gapchars characters that represent gap/missing data (e.g. "-.?N")
 * @param case_sensitive if 0, statechars
 *
 * @return character map
 */
PLL_EXPORT pll_state_t * pllmod_util_charmap_create(unsigned int states,
                                                    const char * statechars,
                                                    const char * gapchars,
                                                    int case_sensitive)
{
  size_t i;
  static const unsigned int maxstates = sizeof(pll_state_t) * 8;

  if (states > maxstates)
  {
    pllmod_set_error(PLLMOD_UTIL_ERROR_MODEL_INVALID_DEF,
                     "The specified number of states (%u) exceeds the allowed maximum (%u)",
                     states, maxstates);
    return NULL;
  }

  if (states > strlen(statechars))
  {
    pllmod_set_error(PLLMOD_UTIL_ERROR_MODEL_INVALID_MAPSTRING,
                     "Character map string is too short for a given number of states: %u",
                     states);
    return NULL;
  }

  pll_state_t * map = calloc(256, sizeof(pll_state_t));

  /* fill map */
  pll_state_t state = 1;
  pll_state_t gapstate = 0;
  for (i = 0; i < states;  ++i)
  {
    int c = statechars[i];

    if (case_sensitive)
    {
      map[c] = state;
    }
    else
    {
      map[tolower(c)] = state;
      map[toupper(c)] = state;
    }

    gapstate |= state;

    state <<= 1;
  }

  assert(( (unsigned int) PLL_STATE_CTZ(state) == states) ||
           (states == maxstates && state == 0));
  assert(PLL_STATE_POPCNT(gapstate) == states);

  /* fill gaps */
  if (gapchars)
  {
    for (i = 0; i < strlen(gapchars);  ++i)
    {
      int c = gapchars[i];
      assert(!map[c]);
      map[c] = gapstate;
    }
  }

  return map;
}

