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

 /* Consensus tree example
  * create a consensus tree out of a set of trees.
  *
  * input: trees_file consensus_treshold
  *        where `trees_file` is a file where each line contains one tree in NEWICK format
  *        for example:
  *        (A,(B,E),(C,D));
  *        ((A,(C,D)),B,E);
  *        `consensus_threshold` [0.0, 1.0] is the minimum branch support
  *
  * output: prints the consensus tree with branch supports.
  */

#include <libpll/pll_tree.h>
#include <assert.h>
#include <stdarg.h>
#if(USE_HASH)
#include <search.h>
#endif
#include <time.h>

static void fatal (const char * format, ...);

/* function for printing non-binary trees */
static void print_newick_recurse(pll_unode_t * node);
static void print_newick(pll_unode_t * tree);

int main (int argc, char * argv[])
{
  unsigned int tree_count;

  if (argc != 3)
    fatal (" syntax: %s [trees_file] [consensus_threshold]", argv[0]);

  /* get arguments */
  char * filename = argv[1];
  double threshold = atof(argv[2]);

  /* build consensus tree */
  pll_consensus_utree_t * constree =
    pllmod_utree_consensus(filename,
                          threshold,
                          &tree_count);
  if (!constree)
    fatal("Error %d: %s\n", pll_errno, pll_errmsg);

  /* print it in NEWICK format */
  print_newick(constree->tree);

  /* clean up */
  pllmod_utree_consensus_destroy(constree);

  return 0;
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

typedef struct consensus_data
{
  pll_split_t split;
  unsigned int bit_count;
  double support;
} consensus_data_t;

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
  if (node->data)
  {
    consensus_data_t * cdata = (consensus_data_t *) node->data;
    printf("[%.3f]", cdata->support);
  }
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
