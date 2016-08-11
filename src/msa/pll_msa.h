/*
 Copyright (C) 2016 Diego Darriba

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
#ifndef PLL_MSA_H_
#define PLL_MSA_H_

#ifndef PLL_H_
#define PLL_H_
#include "pll.h"
#endif

PLL_EXPORT double * pllmod_msa_empirical_frequencies(pll_partition_t * partition);
PLL_EXPORT double * pllmod_msa_empirical_subst_rates(pll_partition_t * partition);
PLL_EXPORT double pllmod_msa_empirical_invariant_sites(pll_partition_t *partition);

#endif /* PLL_MSATREE_H_ */
