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

#include "pll_tree.h"
#include "pll_optimize.h"
#include "pllmod_algorithm.h"
#include <stdarg.h>
#include <search.h>
#include <time.h>

#define STATES    4
#define RATE_CATS 4
#define BRLEN_MIN 1e-6
#define BRLEN_MAX 1e+2

static void fatal(const char * format, ...) __attribute__ ((noreturn));

typedef struct
{
  int clv_valid;
} node_info_t;

static void set_missing_branch_length_recursive(pll_utree_t * tree,
                                                double length)
{
  if (tree)
  {
    /* set branch length to default if not set */
    if (!tree->length)
      tree->length = length;

    if (tree->next)
    {
      if (!tree->next->length)
        tree->next->length = length;

      if (!tree->next->next->length)
        tree->next->next->length = length;

      set_missing_branch_length_recursive(tree->next->back, length);
      set_missing_branch_length_recursive(tree->next->next->back, length);
    }
  }
}

/* branch lengths not present in the newick file get a value of 0.000001 */
static void set_missing_branch_length(pll_utree_t * tree, double length)
{
  set_missing_branch_length_recursive(tree, length);
  set_missing_branch_length_recursive(tree->back, length);
}

static void fatal(const char * format, ...)
{
  va_list argptr;
  va_start(argptr, format);
  vfprintf(stderr, format, argptr);
  va_end(argptr);
  fprintf(stderr, "\n");
  exit(EXIT_FAILURE);
}

int main(int argc, char * argv[])
{
  unsigned int i;
  unsigned int tip_nodes_count, inner_nodes_count, nodes_count, branch_count;
  unsigned int sequence_count;
  pll_partition_t * partition;

  /* we accept only two arguments - the newick tree (unrooted binary) and the
     alignment in the form of FASTA reads */
  if (argc != 3)
    fatal(" syntax: %s [newick] [phylip]", argv[0]);

  /* parse the unrooted binary tree in newick format, and store the number
     of tip nodes in tip_nodes_count */
  pll_utree_t * tree = pll_utree_parse_newick(argv[1], &tip_nodes_count);
  if (!tree)
    fatal("Tree must be an unrooted binary tree");

  /* fix all missing branch lengths (i.e. those that did not appear in the
     newick) to 0.000001 */
  set_missing_branch_length(tree, BRLEN_MIN);

  /* compute and show node count information */
  inner_nodes_count = tip_nodes_count - 2;
  nodes_count = inner_nodes_count + tip_nodes_count;
  branch_count = nodes_count - 1;

  printf("Number of tip/leaf nodes in tree: %d\n", tip_nodes_count);
  printf("Number of inner nodes in tree: %d\n", inner_nodes_count);
  printf("Total number of nodes in tree: %d\n", nodes_count);
  printf("Number of branches in tree: %d\n\n", branch_count);

  /* Uncomment to display the parsed tree ASCII tree together with information
     as to which CLV index, branch length and label is associated with each
     node. The code will also write (and print on screen) the newick format
     of the tree.

  pll_utree_show_ascii(tree, PLL_UTREE_SHOW_LABEL |
                             PLL_UTREE_SHOW_BRANCH_LENGTH |
                             PLL_UTREE_SHOW_CLV_INDEX);
  char * newick = pll_utree_export_newick(tree);
  printf("%s\n", newick);
  free(newick);

  */

  /*  obtain an array of pointers to tip nodes */
  pll_utree_t ** tipnodes = (pll_utree_t  **)calloc(tip_nodes_count,
                                                    sizeof(pll_utree_t *));
  pll_utree_query_tipnodes(tree, tipnodes);

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
  pll_msa_t * msa = pll_phylip_parse_msa(argv[2], &sequence_count);
  if (!msa)
    fatal(pll_errmsg);

  /* compress site patterns */
  if (sequence_count != tip_nodes_count)
    fatal("Number of sequences does not match number of leaves in tree");

  printf("Original sequence (alignment) length : %d\n", msa->length);
  unsigned int * weight = pll_compress_site_patterns(msa->sequence,
                                                     pll_map_nt,
                                                     tip_nodes_count,
                                                     &(msa->length));
  printf("Number of unique site patterns: %d\n\n", msa->length);


  /* create the PLL partition instance

  tip_nodes_count : the number of tip sequences we want to have
  inner_nodes_count : the number of CLV buffers to be allocated for inner nodes
  STATES : the number of states that our data have
  1 : number of different substitution models (or eigen decomposition)
      to use concurrently (i.e. 4 for LG4)
  branch_count: number of probability matrices to be allocated
  RATE_CATS : number of rate categories we will use
  inner_nodes_count : how many scale buffers to use
  PLL_ATTRIB_ARCH_SSE : list of flags for hardware acceleration (not yet implemented)

  */

  partition = pll_partition_create(tip_nodes_count,
                                   inner_nodes_count,
                                   STATES,
                                   (unsigned int)(msa->length),
                                   1,
                                   branch_count,
                                   RATE_CATS,
                                   inner_nodes_count,
                                   PLL_ATTRIB_ARCH_AVX);

  /* initialize the array of base frequencies */
  double frequencies[4] = { 0.17, 0.19, 0.25, 0.39 };

  /* substitution rates for the 4x4 GTR model. This means we need exactly
     (4*4-4)/2 = 6 values, i.e. the number of elements above the diagonal */
  double subst_params[6] = {1,1,1,1,1,1};

  /* we'll use 4 rate categories, and currently initialize them to 0 */
  double rate_cats[4] = {0};

  /* compute the discretized category rates from a gamma distribution
     with alpha shape 1 and store them in rate_cats  */
  pll_compute_gamma_cats(1, 4, rate_cats);

  /* set frequencies at model with index 0 (we currently have only one model) */
  pll_set_frequencies(partition, 0, frequencies);

  /* set 6 substitution parameters at model with index 0 */
  pll_set_subst_params(partition, 0, subst_params);

  /* set rate categories */
  pll_set_category_rates(partition, rate_cats);

  /* set pattern weights and free the weights array */
  pll_set_pattern_weights(partition, weight);
  free(weight);

  /* find sequences in hash table and link them with the corresponding taxa */
  for (i = 0; i < tip_nodes_count; ++i)
  {
    ENTRY query;
    query.key = msa->label[i];
    ENTRY * found = NULL;

    found = hsearch(query,FIND);

    if (!found)
      fatal("Sequence with header %s does not appear in the tree", msa->label[i]);

    unsigned int tip_clv_index = *((unsigned int *)(found->data));

    pll_set_tip_states(partition, tip_clv_index, pll_map_nt, msa->sequence[i]);
  }

  pll_msa_destroy(msa);


  /* destroy hash table */
  hdestroy();

  /* we no longer need these two arrays (keys and values of hash table... */
  free(data);
  free(tipnodes);


  /* update matrix_count probability matrices using the rate matrix with
     index 0. The i-th matrix (i ranges from 0 to matrix_count - 1) is
     generated using branch length branch_lengths[i] and rate matrix
     (substitution rates + frequencies) params_indices[i], and can be refered
     to with index matrix_indices[i] */
  unsigned int params_indices[4] = {0,0,0,0};

  int params_to_optimize = PLLMOD_OPT_PARAM_BRANCHES_ITERATIVE;

  /* create treeinfo structure */
  pllmod_treeinfo_t * treeinfo = pllmod_treeinfo_create(tree,
                                                        tip_nodes_count,
                                                        1,
                                                        PLLMOD_TREE_BRLEN_LINKED);

  int retval = pllmod_treeinfo_init_partition(treeinfo,
                                              0,
                                              partition,
                                              params_to_optimize,
                                              1.0, /* alpha*/
                                              params_indices, /* param_indices */
                                              NULL /* subst matrix symmetries*/
                                              );

  if (!retval)
    fatal("Error initializing partition!");

  /* Compute initial LH of the starting tree */
  double loglh = pllmod_treeinfo_compute_loglh(treeinfo, 0);

  printf("Log-L before SPRs: %lf\n", loglh);

  /* define SPR parameters */
  int spr_thorough = 0;          /* perform triplet BLO after each SPR? */
  int spr_radius_min = 1;        /* MIN re-insertion radius */
  int spr_radius_max = 5;        /* MAX re-insertion radius */
  int spr_ntopol_keep = 20;      /* topologies for SLOW re-evaluation (full BLO) */
  double spr_subtree_cutoff = 0.; /* not used here */
  double spr_lh_epsilon = 0.1;   /* logLH epsilon */
  double spr_blo_smoothings = 8; /* MAX number of BLO iterations */

  /* no do a round of SPRs */
  loglh = pllmod_algo_spr_round(treeinfo,
                                spr_radius_min,
                                spr_radius_max,
                                spr_ntopol_keep,
                                spr_thorough,
                                BRLEN_MIN,
                                BRLEN_MAX,
                                spr_blo_smoothings,
                                spr_lh_epsilon,
                                NULL,                /* cutoff_info: not used here */
                                spr_subtree_cutoff
                               );

  printf("Log-L after SPRs: %lf\n\n", loglh);

  pllmod_treeinfo_destroy(treeinfo);

  /* destroy all structures allocated for the concrete PLL partition instance */
  pll_partition_destroy(partition);

  /* we will no longer need the tree structure */
  pll_utree_destroy(treeinfo->root);

  return (EXIT_SUCCESS);
}
