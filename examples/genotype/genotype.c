/*
    Copyright (C) 2019 Alexey Kozlov

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

#include "libpll/pll_tree.h"
#include "libpll/pll_optimize.h"
#include "libpll/pllmod_algorithm.h"
#include "libpll/pllmod_common.h"
#include <stdarg.h>
#include <search.h>
#include <time.h>

#define GT_MODEL "GTGTR4"
//#define GT_MODEL "GTJC"
#define RATE_CATS 1
#define BRLEN_MIN 1e-6
#define BRLEN_MAX 1e+2

static void fatal(const char * format, ...) __attribute__ ((noreturn));

static void fatal(const char * format, ...)
{
  va_list argptr;
  va_start(argptr, format);
  vfprintf(stderr, format, argptr);
  va_end(argptr);
  fprintf(stderr, "\n");
  exit(EXIT_FAILURE);
}

pll_utree_t * load_tree(const char *fname)
{
  /* parse the unrooted binary tree in newick format, and store the number
     of tip nodes in tip_nodes_count */
  pll_utree_t * tree = pll_utree_parse_newick(fname);
  if (!tree)
    fatal("Tree must be an unrooted binary tree");

  /* fix all missing branch lengths (i.e. those that did not appear in the
     newick) to 0.000001 */
  pllmod_utree_set_length_recursive(tree, BRLEN_MIN, 1);

  return tree;
}

void set_partition_tips(pll_partition_t * partition, pll_msa_t * msa)
{
  int i;

  /* find sequences in hash table and link them with the corresponding taxa */
  for (i = 0; i < msa->count; ++i)
  {
    ENTRY query;
    query.key = msa->label[i];
    ENTRY * found = NULL;

    found = hsearch(query,FIND);

    if (!found)
      fatal("Sequence with header %s does not appear in the tree", msa->label[i]);

    unsigned int tip_clv_index = *((unsigned int *)(found->data));

    pll_set_tip_states(partition, tip_clv_index, pll_map_gt10, msa->sequence[i]);
  }
}

double * expand_uniq_rates(int states, const double * uniq_rates, const int * rate_sym)
{
  unsigned int i;

  unsigned int num_rates = states * (states-1) / 2;
  double * subst_rates = calloc(num_rates, sizeof(double));
  for (i = 0; i < num_rates; ++i)
    subst_rates[i] = rate_sym ? uniq_rates[rate_sym[i]] : uniq_rates[i];

  return subst_rates;
}

int main(int argc, char * argv[])
{
  unsigned int i;
  unsigned int tip_nodes_count, inner_nodes_count, nodes_count, branch_count;
  pll_partition_t * partition;

  /* we accept only two arguments - the newick tree (unrooted binary) and the
     alignment in the form of FASTA reads */
  if (argc != 3)
    fatal(" syntax: %s [newick] [phylip]", argv[0]);

  pll_utree_t * tree = load_tree(argv[1]);

  /* compute and show node count information */
  tip_nodes_count = tree->tip_count;
  inner_nodes_count = tree->inner_count;
  nodes_count = inner_nodes_count + tip_nodes_count;
  branch_count = tree->edge_count;

  printf("Number of tip/leaf nodes in tree: %d\n", tip_nodes_count);
  printf("Number of inner nodes in tree: %d\n", inner_nodes_count);
  printf("Total number of nodes in tree: %d\n", nodes_count);
  printf("Number of branches in tree: %d\n\n", branch_count);

  /* obtain an array of pointers to tip nodes:
   * they are always stored at the beginning of the tree->nodes array  */
  pll_unode_t ** tipnodes = tree->nodes;

  /* create a libc hash table of size tip_nodes_count */
  hcreate(tip_nodes_count);

  /* populate a libc hash table with tree tip labels */
  unsigned int * data = (unsigned int *)malloc(tip_nodes_count *
                                               sizeof(unsigned int));
  for (i = 0; i < tip_nodes_count; ++i)
  {
    data[i] = tipnodes[i]->clv_index;
    ENTRY entry;
    entry.key = tipnodes[i]->label;
    entry.data = (void *)(data+i);
    hsearch(entry, ENTER);
  }

  /* read PHYLIP alignment */
  pll_msa_t * msa = pll_phylip_load(argv[2], PLL_FALSE);
  if (!msa)
    fatal(pll_errmsg);

  /* compress site patterns */
  if (msa->count != (int) tip_nodes_count)
    fatal("Number of sequences does not match number of leaves in tree");

  printf("Original sequence (alignment) length : %d\n", msa->length);
  unsigned int * weight = pll_compress_site_patterns(msa->sequence,
                                                     pll_map_gt10,
                                                     tip_nodes_count,
                                                     &(msa->length));
  printf("Number of unique site patterns: %d\n\n", msa->length);



  pllmod_subst_model_t * model = pllmod_util_model_info_genotype(GT_MODEL);

  if (!model)
    fatal("Unknown evolutionary model: %s", GT_MODEL);

  /* create the PLL partition instance

  tip_nodes_count : the number of tip sequences we want to have
  inner_nodes_count : the number of CLV buffers to be allocated for inner nodes
  model->states : the number of states that our data have
  1 : number of different substitution models (or eigen decomposition)
      to use concurrently (i.e. 4 for LG4)
  branch_count: number of probability matrices to be allocated
  RATE_CATS : number of rate categories we will use
  inner_nodes_count : how many scale buffers to use
  PLL_ATTRIB_ARCH_AVX : list of flags for hardware acceleration

  */

  partition = pll_partition_create(tip_nodes_count,
                                   inner_nodes_count,
                                   model->states,
                                   (unsigned int)(msa->length),
                                   1,
                                   branch_count,
                                   RATE_CATS,
                                   inner_nodes_count,
                                   PLL_ATTRIB_ARCH_AVX);


  /* initialize the array of base frequencies */
  double user_freqs[10] = { 0.254504, 0.157238, 0.073667, 0.287452, 0.027336,
                            0.058670, 0.041628, 0.024572, 0.064522, 0.010412 };

  /* substitution rates: for GTR4 model those are 6 "regular" DNA susbt. rates + 1 rate
   * for "unlikely" double substitutions (eg A/A -> C/T) */
  double unique_subst_rates[7] = { 0.001000, 0.101223, 0.001000, 0.001000, 1.000000,
                                   0.001000, 0.447050 };

  /* get full above-diagonal half-matrix */
  double * user_subst_rates = expand_uniq_rates(model->states, unique_subst_rates,
                                                model->rate_sym);

  double rate_cats[RATE_CATS] = {0};

  /* compute the discretized category rates from a gamma distribution
     with alpha shape 1 and store them in rate_cats  */
  pll_compute_gamma_cats(1, RATE_CATS, rate_cats, PLL_GAMMA_RATES_MEAN);

  /* set frequencies at model with index 0 (we currently have only one model) */
  pll_set_frequencies(partition, 0, model->freqs ? model->freqs : user_freqs);

  /* set substitution parameters at model with index 0 */
  pll_set_subst_params(partition, 0, model->rates ? model->rates : user_subst_rates);
  free(user_subst_rates);

  /* set rate categories */
  pll_set_category_rates(partition, rate_cats);

  /* set pattern weights and free the weights array */
  pll_set_pattern_weights(partition, weight);
  free(weight);

  set_partition_tips(partition, msa);

  pll_msa_destroy(msa);

  /* destroy hash table */
  hdestroy();

  /* we no longer need these two arrays (keys and values of hash table... */
  free(data);


  /* update matrix_count probability matrices using the rate matrix with
     index 0. The i-th matrix (i ranges from 0 to matrix_count - 1) is
     generated using branch length branch_lengths[i] and rate matrix
     (substitution rates + frequencies) params_indices[i], and can be refered
     to with index matrix_indices[i] */
  unsigned int params_indices[RATE_CATS] = {0};

  /* we do not want to optimize anything */
  int params_to_optimize = 0;

  /* create treeinfo structure */
  pllmod_treeinfo_t * treeinfo = pllmod_treeinfo_create(tree->vroot,
                                                        tip_nodes_count,
                                                        1,
                                                        PLLMOD_COMMON_BRLEN_LINKED);

  int retval = pllmod_treeinfo_init_partition(treeinfo,
                                              0,
                                              partition,
                                              params_to_optimize,
                                              PLL_GAMMA_RATES_MEAN,
                                              1.0, /* alpha*/
                                              params_indices, /* param_indices */
                                              model->rate_sym /* subst matrix symmetries*/
                                              );


  if (!retval)
    fatal("Error initializing partition!");

  /* Compute initial LH of the starting tree */
  double loglh = pllmod_treeinfo_compute_loglh(treeinfo, 0);

  printf("Log-Likelihood : %lf\n", loglh);

  pllmod_treeinfo_destroy(treeinfo);

  /* destroy all structures allocated for the concrete PLL partition instance */
  pll_partition_destroy(partition);

  /* we will no longer need the tree structure */
  pll_utree_destroy(tree, NULL);

  pllmod_util_model_destroy(model);

  return (EXIT_SUCCESS);
}
