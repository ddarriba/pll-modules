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
#include "pll_tree.h"
#include "../rng.h"
#include "../common.h"

#include <assert.h>

#define TREEFILE  "testdata/medium.tree"

/*
 * This example reads an input tree, generates the set of splits, shuffles them
 * randomly and reconstructs the tree out of the splits
 */
/* static functions */
static void shuffle(pll_split_t *array, size_t n);
static void print_newick_recurse(pll_unode_t * node);
static void print_newick(pll_unode_t * tree);

static const unsigned int n_iters = 10;
int main (int argc, char * argv[])
{
  unsigned int rf_dist, i, iter;
  char **labels;

  /* tree properties */
  pll_unode_t * tree = NULL;
  unsigned int tip_count;
  unsigned int attributes = get_attributes(argc, argv);

  if (attributes != PLL_ATTRIB_ARCH_CPU)
  {
    skip_test();
  }

  /* parse the input trees */
  pll_utree_t * parsed_tree = pll_utree_parse_newick (TREEFILE);
  tip_count = parsed_tree->tip_count;
  tree = parsed_tree->nodes[2*tip_count - 3];
  if (!tree)
  {
    fatal("Error %d: %s", pll_errno, pll_errmsg);
  }

  pll_unode_t ** tipnodes = parsed_tree->nodes;

  labels = (char **) malloc(tip_count * sizeof(char *));
  for (i=0; i<tip_count; ++i)
    labels[tipnodes[i]->node_index] = tipnodes[i]->label;

  unsigned int n_splits = tip_count - 3;
  pll_split_t * splits = pllmod_utree_split_create(tree,
                                                   tip_count,
                                                   NULL);

  for (iter=0; iter<n_iters; ++iter)
  {
    shuffle(splits, n_splits);

    pll_split_system_t split_system;
    split_system.splits = splits;
    split_system.support = 0;
    split_system.split_count = n_splits;
    split_system.max_support = 1.0;

    pll_consensus_utree_t * constree = pllmod_utree_from_splits(&split_system,
                                                                tip_count,
                                                                labels);

    pll_utree_t * consensus = pll_utree_wraptree(constree->tree, tip_count);
    if (!pllmod_utree_consistency_set(consensus, parsed_tree))
       fatal("Cannot set trees consistent!");

    pll_utree_show_ascii(constree->tree, PLL_UTREE_SHOW_CLV_INDEX | PLL_UTREE_SHOW_LABEL);
    print_newick(constree->tree);

    pll_split_t * splits2 = pllmod_utree_split_create(constree->tree,
                                                      tip_count,
                                                      NULL);

    pllmod_utree_split_normalize_and_sort(splits2,
                                          tip_count,
                                          n_splits,
                                          1);

    /* sort splits back */
    pllmod_utree_split_normalize_and_sort(splits,
                                          tip_count,
                                          n_splits,
                                          0);

    rf_dist = pllmod_utree_split_rf_distance(splits, splits2, tip_count);
    printf(" RF DIST = %d\n", rf_dist);

    if (rf_dist > 0)
      fatal("Error: Initial and reconstructed trees differ!");

    /* in-loop cleanup */
    pllmod_utree_consensus_destroy(constree);
    free (consensus->nodes);
    free (consensus);
    pllmod_utree_split_destroy(splits2);
  }

  /* clean */
  free(labels);
  pll_utree_destroy (parsed_tree, NULL);
  pllmod_utree_split_destroy(splits);

  return (0);
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

static void shuffle(pll_split_t *array, size_t n)
{
  unsigned int i, j;
  pll_split_t t;
  for (i = 0; i < n - 1; i++)
  {
    j = i + (unsigned int) RAND / (RAND_MAX / (n - i) + 1);
    t = array[j];
    array[j] = array[i];
    array[i] = t;
  }
}

static void print_newick_recurse(pll_unode_t * node)
{
  pll_unode_t * child;
  if (pllmod_utree_is_tip(node))
  {
    printf("%s", node->label);
    return;
  }

  printf("(");
  child = node->next;
  while(child != node)
  {
    print_newick_recurse(child->back);

    if (child->next != node)
      printf(",");

    child = child->next;
  }
  printf(")");
}

static void print_newick(pll_unode_t * tree)
{
  printf("(");
  print_newick_recurse(tree->back);
  pll_unode_t * child = tree->next;
  while(child != tree)
  {
    printf(",");
    print_newick_recurse(child->back);
    child = child->next;
  }
  printf(");\n");
}
