/*
 Copyright (C) 2016 Diego Darriba, Alexey Kozlov

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

#ifndef PLL_MOD_UTIL_H_
#define PLL_MOD_UTIL_H_

#ifndef PLL_H_
#define PLL_H_
#include "pll.h"
#endif

/* number of substitution rates in protein model matrices */
#define PLLMOD_PROT_RATES_COUNT 190

/* number of aminoacid frequencies in model specification */
#define PLLMOD_PROT_FREQS_COUNT 20

/* error codes for UTIL libpll module (5001-6000)*/
#define PLLMOD_ERROR_UNKNOWN_MODEL                5001


/* functions for working with built-in protein models */
PLL_EXPORT unsigned int pllmod_util_get_protein_model_count();
PLL_EXPORT const char ** pllmod_util_get_protein_model_names();
PLL_EXPORT int pllmod_util_protein_model_exists(const char * model_name);
PLL_EXPORT const double * pllmod_util_get_protein_model_rates(const char * model_name);
PLL_EXPORT const double * pllmod_util_get_protein_model_freqs(const char * model_name);
PLL_EXPORT unsigned int pllmod_util_get_protein_model_matrix_count(const char * model_name);
PLL_EXPORT int pllmod_util_set_protein_model(pll_partition_t * partition, const char * model_name, int model_freqs);

#endif
