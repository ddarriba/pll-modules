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
#include <libpll/pll_tree.h>

#include <assert.h>
#include <stdarg.h>

/* set to 1 for printing splits */
#define PRINT_SPLITS 0

/* static functions */
static void fatal (const char * format, ...);

int main (int argc, char * argv[])
{
  unsigned int rf_dist;

  /* tree properties */
  pll_utree_t * tree1 = NULL,
              * tree2 = NULL;
  unsigned int tip_count;

  if (argc != 3)
    fatal (" syntax: %s [newick] [newick]", argv[0]);

  /* parse the input trees */
  tree1 = pll_utree_parse_newick (argv[1]);
  tree2 = pll_utree_parse_newick (argv[2]);
  tip_count = tree1->tip_count;

  if (tip_count != tree2->tip_count)
    fatal("Trees have different number of tips!");

  if (!pllmod_utree_consistency_set(tree1, tree2))
    fatal("Cannot set trees consistent!");

  if (!pllmod_utree_consistency_check(tree1, tree2))
    fatal("Tip node IDs are not consistent!");

  /* uncomment lines below for displaying the trees in ASCII format */
  // pll_utree_show_ascii(tree1, PLL_UTREE_SHOW_LABEL | PLL_UTREE_SHOW_PMATRIX_INDEX);
  // pll_utree_show_ascii(tree2, PLL_UTREE_SHOW_LABEL | PLL_UTREE_SHOW_PMATRIX_INDEX);

  /* next, we compute the RF distance in 2 different ways: */

  /* 1. creating the split sets manually */
  unsigned int n_splits = tip_count - 3;
  pll_split_t * splits1 = pllmod_utree_split_create(tree1->nodes[tip_count],
                                                    tip_count,
                                                    NULL);

#if(PRINT_SPLITS)
  {
    unsigned int i;
    for (i=0; i<n_splits; ++i)
    {
      pllmod_utree_split_show(splits1[i], tip_count);
      printf("\n");
    }
    printf("\n");
  }
#endif

  /* now we compute the splits, but also the nodes corresponding to each split */
  pll_unode_t ** splits_to_node = (pll_unode_t **) malloc(n_splits * sizeof(pll_unode_t *));
  pll_split_t * splits2 = pllmod_utree_split_create(tree2->nodes[tip_count],
                                                    tip_count,
                                                    splits_to_node);

#if(PRINT_SPLITS)
  {
    unsigned int i;
    for (i=0; i<n_splits; ++i)
    {
      pllmod_utree_split_show(splits2[i], tip_count);
      printf(" node: Pmatrix:%d Nodes:%d<->%d\n",
                                splits_to_node[i]->pmatrix_index,
                                splits_to_node[i]->node_index,
                                splits_to_node[i]->back->node_index);
    }
  }
#endif

  rf_dist = pllmod_utree_split_rf_distance(splits1, splits2, tip_count);
  printf("RF [manual]\n");
  printf("distance = %d\n", rf_dist);
  printf("relative = %.2f%%\n", 100.0*rf_dist/(2*(tip_count-3)));

  pllmod_utree_split_destroy(splits1);
  pllmod_utree_split_destroy(splits2);
  free(splits_to_node);

  /* 2. directly from the tree structures */

  rf_dist = pllmod_utree_rf_distance(tree1->nodes[tip_count],
                                     tree2->nodes[tip_count],
                                     tip_count);

  printf("RF [auto]\n");
  printf("distance = %d\n", rf_dist);
  printf("relative = %.2f%%\n", 100.0*rf_dist/(2*(tip_count-3)));

  /* clean */
  pll_utree_destroy (tree1, NULL);
  pll_utree_destroy (tree2, NULL);

  return (0);
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

static void fatal (const char * format, ...)
{
  va_list argptr;
  va_start(argptr, format);
  vfprintf (stderr, format, argptr);
  va_end(argptr);
  fprintf (stderr, "\n");
  exit (EXIT_FAILURE);
}
