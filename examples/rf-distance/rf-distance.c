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
//#include "utils.h"
#include "pll_tree.h"
#include <assert.h>
#include <stdarg.h>

/* static functions */
static void fatal (const char * format, ...);
static int cb_set_node_ids(pll_utree_t * node, void *data);
static int cb_clean_data(pll_utree_t * node, void *data);

typedef struct
{
  unsigned int node_id;
  unsigned int subnode_id;
} t_tip_info;

int main (int argc, char * argv[])
{
  unsigned int i,
               rf_dist;

  /* tree properties */
  pll_utree_t * tree1 = NULL,
              * tree2 = NULL,
              ** tipnodes,
              ** tipnodes_alt;
  unsigned int tip_count,
               tip_count_alt,
               nodes_count;

  if (argc != 3)
    fatal (" syntax: %s [newick] [newick]", argv[0]);

  /* parse the input trees */
  tree1 = pll_utree_parse_newick (argv[1], &tip_count);
  tree2 = pll_utree_parse_newick (argv[2], &tip_count_alt);

  if (tip_count != tip_count_alt)
    fatal("Trees have different number of tips!");

  nodes_count = 2 * tip_count - 2;

  /* IMPORTANT NOTE!
   * While the node and subnode ids are not set in pll, we need to manually
   * allocate space for them in the utree data pointer and set them!!
   *
   * These lines below can be removed once this is implemented.
   */

  /*  obtain an array of pointers to tip nodes */

  tipnodes = (pll_utree_t **) calloc ((size_t) tip_count,
                                      sizeof(pll_utree_t *));
  pll_utree_query_tipnodes (tree1, tipnodes);
  tipnodes_alt = (pll_utree_t **) calloc ((size_t) tip_count,
                                          sizeof(pll_utree_t *));
  pll_utree_query_tipnodes (tree2, tipnodes_alt);

  for (i = 0; i < tip_count; ++i)
  {
    tipnodes[i]->data = (t_tip_info*) malloc(sizeof(t_tip_info));
    tipnodes_alt[i]->data = (t_tip_info*) malloc(sizeof(t_tip_info));
    ((t_tip_info *)tipnodes[i]->data)->node_id = i;
  }

  /* set inner node ids in postorder traversal */

  t_tip_info tip_info;
  tip_info.node_id = tip_count;
  pll_utree_traverse_apply(tree1,
                           0,
                           &cb_set_node_ids,
                           &tip_info);
  if(tip_info.node_id != nodes_count)
    fatal("Number of nodes do not agree!");
  tip_info.node_id = tip_count;
  pll_utree_traverse_apply(tree2,
                           0,
                           &cb_set_node_ids,
                           &tip_info);
  if(tip_info.node_id != nodes_count)
    fatal("Number of nodes do not agree!");

  if (!pll_utree_consistency_set(tree1, tree2, tip_count))
    fatal("Cannot set trees consistent!");

  if (!pll_utree_consistency_check(tree1, tree2, tip_count))
    fatal("Tip node IDs are not consistent!");

  /* uncomment lines below for displaying the trees in ASCII format */
  // pll_utree_show_ascii(tree1, PLL_UTREE_SHOW_LABEL | PLL_UTREE_SHOW_PMATRIX_INDEX);
  // pll_utree_show_ascii(tree2, PLL_UTREE_SHOW_LABEL | PLL_UTREE_SHOW_PMATRIX_INDEX);

  /* next, we compute the RF distance in 2 different ways: */

  /* 1. creating the split sets manually */
  unsigned int n_splits;
  pll_split_t * splits1 = pll_utree_split_create(tree1,
                                                 tip_count,
                                                 &n_splits);

  pll_split_t * splits2 = pll_utree_split_create(tree2,
                                                 tip_count,
                                                 &n_splits);

  /* uncomment lines below for printing the splits out */
  // printf("\nSPLITS 1\n");
  // for (i=0; i<n_splits; ++i)
  //   pll_utree_split_show(splits1[i], tip_count);
  // printf("\nSPLITS 2\n");
  // for (i=0; i<n_splits; ++i)
  //     pll_utree_split_show(splits2[i], tip_count);
  // printf("\n");

  rf_dist = pll_utree_split_rf_distance(splits1, splits2, tip_count);
  printf("RF [manual]\n");
  printf("distance = %d\n", rf_dist);
  printf("relative = %.2f%%\n", 100.0*rf_dist/(2*(tip_count-3)));

  pll_utree_split_destroy(splits1);
  pll_utree_split_destroy(splits2);

  /* 2. directly from the tree structures */

  rf_dist = pll_utree_rf_distance(tree1,
                                  tree2,
                                  tip_count);

  printf("RF [auto]\n");
  printf("distance = %d\n", rf_dist);
  printf("relative = %.2f%%\n", 100.0*rf_dist/(2*(tip_count-3)));

  /* clean */
  pll_utree_traverse_apply(tree1, 0, &cb_clean_data, 0);
  pll_utree_destroy (tree1);

  pll_utree_traverse_apply(tree2, 0, &cb_clean_data, 0);
  pll_utree_destroy (tree2);

  free (tipnodes);
  free (tipnodes_alt);

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

static int cb_set_node_ids(pll_utree_t * node, void *data)
{
  /* assign and increase the node id, assign subnode ids from 0 to 3 */
  if (node->next)
  {
    t_tip_info * tip_info = (t_tip_info*) data;
    node->data =  malloc(sizeof(t_tip_info));
    ((t_tip_info *)node->data)->node_id = tip_info->node_id;
    ((t_tip_info *)node->data)->subnode_id = 0;

    node->next->data =  malloc(sizeof(t_tip_info));
    ((t_tip_info *)node->next->data)->node_id = tip_info->node_id;
    ((t_tip_info *)node->next->data)->subnode_id = 1;

    node->next->next->data =  malloc(sizeof(t_tip_info));
    ((t_tip_info *)node->next->next->data)->node_id = tip_info->node_id;
    ((t_tip_info *)node->next->next->data)->subnode_id = 2;

    ++tip_info->node_id;
  }

  /* continue traversing */
  return 1;
}

static int cb_clean_data(pll_utree_t * node, void *data)
{
  /* clean the void data pointers */
  free(node->data);
  if (node->next)
  {
    free(node->next->data);
    free(node->next->next->data);
  }

  /* continue traversing */
  return 1;
}
