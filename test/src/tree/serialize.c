/*
 Copyright (C) 2015 Diego Darriba, Tomas Flouri

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
#include "../common.h"

#include <string.h>
#include <stdarg.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>

#define N_TAXA_SMALL 100

const char * header[100] = {
  "T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10",
  "T11", "T12", "T13", "T14", "T15", "T16", "T17", "T18", "T19", "T20",
  "T21", "T22", "T23", "T24", "T25", "T26", "T27", "T28", "T29", "T30",
  "T31", "T32", "T33", "T34", "T35", "T36", "T37", "T38", "T39", "T40",
  "T41", "T42", "T43", "T44", "T45", "T46", "T47", "T48", "T49", "T50",
  "T51", "T52", "T53", "T54", "T55", "T56", "T57", "T58", "T59", "T60",
  "T61", "T62", "T63", "T64", "T65", "T66", "T67", "T68", "T69", "T70",
  "T71", "T72", "T73", "T74", "T75", "T76", "T77", "T78", "T79", "T80",
  "T81", "T82", "T83", "T84", "T85", "T86", "T87", "T88", "T89", "T90",
  "T91", "T92", "T93", "T94", "T95", "T96", "T97", "T98", "T99", "T100"};

static int cb_set_bl(pll_unode_t * tree,
                     void * data)
{
  assert(tree->pmatrix_index == tree->back->pmatrix_index);
  tree->length = tree->back->length = 1.0*tree->pmatrix_index;
  return 1;
}

static int cb_set_names(pll_unode_t * tree,
                        void * data)
{
  const char ** names = (const char **) data;
  int tip_index;

  if (!tree->next)
  {
    tip_index = tree->node_index;
    tree->label = (char *) malloc (strlen(names[tip_index]) + 1);
    strcpy(tree->label, names[tip_index]);
  }
  return 1;
}

static pll_unode_t * get_utree_root(pll_utree_t * tree)
{
  return tree->nodes[tree->tip_count + tree->inner_count - 1];
}

int main (int argc, char * argv[])
{
   unsigned int n_taxa = N_TAXA_SMALL;
   //unsigned int attributes = get_attributes(argc, argv);

   /* fix RNG seed to arbitrary number */
   srand(42);

   pll_utree_t * random_tree = pllmod_utree_create_random(n_taxa,
                                                         (const char **)header);
   pll_unode_t * root = random_tree->nodes[2*n_taxa - 3];

   /* set arbitrary branch lengths */
   pllmod_utree_traverse_apply(root,
                               NULL,
                               NULL,
                               cb_set_bl,
                               NULL);

   printf("Root set to %d\n", root->node_index);

   printf("\nINITIAL RANDOM TREE:\n\n");
   pll_utree_show_ascii(root, PLL_UTREE_SHOW_CLV_INDEX | PLL_UTREE_SHOW_BRANCH_LENGTH | PLL_UTREE_SHOW_LABEL | PLL_UTREE_SHOW_PMATRIX_INDEX);

   int rf_distance;
   pll_utree_t * tree;
   pll_utree_t * tree2;
   pll_unode_t * stree;
   pll_unode_t * root2;

   /* serialize and expand */
   stree = pllmod_utree_serialize(root, n_taxa);
   tree2 = pllmod_utree_expand(stree, n_taxa);
   root2 = get_utree_root(tree2);
   free(stree);

   /* reset names */
   pllmod_utree_traverse_apply(root2,
                               NULL,
                               NULL,
                               cb_set_names,
                               header);

   printf("\n\nRECONSTRUCTED TREE:\n\n");
   pll_utree_show_ascii(root2, PLL_UTREE_SHOW_CLV_INDEX | PLL_UTREE_SHOW_BRANCH_LENGTH | PLL_UTREE_SHOW_LABEL | PLL_UTREE_SHOW_PMATRIX_INDEX);

   rf_distance = pllmod_utree_rf_distance(root,
                                          root2,
                                          n_taxa);
   assert(!rf_distance);

   /* try now with root at tip */
   //pll_utree_destroy(tree, NULL);
   while(root2->node_index > n_taxa) root2 = root2->next?root2->next->back:root2->back;
   printf("Root set to %d\n", root2->node_index);

   stree = pllmod_utree_serialize(root2->back, n_taxa);
   tree = pllmod_utree_expand(stree, n_taxa);
   root = get_utree_root(tree);
   free(stree);

   /* reset names */
   pllmod_utree_traverse_apply(root,
                               NULL,
                               NULL,
                               cb_set_names,
                               header);

   printf("\nRECONSTRUCTED FROM TIP:\n\n");
   pll_utree_show_ascii(root, PLL_UTREE_SHOW_CLV_INDEX | PLL_UTREE_SHOW_BRANCH_LENGTH | PLL_UTREE_SHOW_LABEL | PLL_UTREE_SHOW_PMATRIX_INDEX);

   rf_distance = pllmod_utree_rf_distance(root,
                                          root2,
                                          n_taxa);
   assert(!rf_distance);

   pll_utree_destroy(random_tree, NULL);
   pll_utree_destroy(tree, NULL);
   pll_utree_destroy(tree2, NULL);

   return PLL_SUCCESS;
}
