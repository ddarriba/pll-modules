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

/* Weighted consensus tree example
 * create a consensus tree out of a weighted set of trees.
 *
 * input: trees_file consensus_treshold
 *        `trees_file` is a file where each line has the following format:
 *
 *        weight newick_tree;
 *
 *        for example:
 *        0.4 (A,(B,E),(C,D));
 *        0.6 ((A,(C,D)),B,E);
 *
 *        `consensus_threshold` [0.0, 1.0] is the minimum branch support
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

static FILE *get_number_of_trees(unsigned int *tree_count,
                                 const char *filename)
{
  FILE
    *f = fopen(filename, "r");

  if (!f)
    return NULL;

  unsigned int trees = 0;
  int ch;

  while((ch = fgetc(f)) != EOF)
    if(ch == ';')
      trees++;

  *tree_count = trees;

  rewind(f);

  return f;
}

int main (int argc, char * argv[])
{
  unsigned int tree_count;
  pll_utree_t ** trees; /* set of tree structures */
  double * weights;     /* weight (or support) for each tree */

  if (argc != 3)
    fatal (" syntax: %s [trees_file] [consensus_threshold]", argv[0]);

  /* parse arguments */
  char * filename = argv[1];
  double threshold = atof(argv[2]);

  /* open file and count number of entries */
  FILE *f = get_number_of_trees(&tree_count,
                                filename);
  if (!f)
    fatal("Cannot open %s for reading\n", filename);
  if (!tree_count)
    fatal("File %s does not contain any tree\n", filename);

  trees = (pll_utree_t **) malloc(tree_count * sizeof(pll_utree_t *));
  weights = (double *) malloc(tree_count * sizeof(double));

  /* parse trees of up to 1,000 characters */
  char tree[1000];
  unsigned int cur_tree = 0;
  while(fscanf(f, "%lf %s\n", weights+cur_tree, tree) != -1)
  {
    trees[cur_tree] = pll_utree_parse_newick_string(tree);
    if (pll_errno)
      fatal("Error: %s\n", pll_errmsg);

    /* set node indices consistent with each other */
    if (cur_tree)
      pllmod_utree_consistency_set(trees[0], trees[cur_tree]);
    ++cur_tree;
  }
  fclose(f);
  assert(cur_tree == tree_count);

  /* build consensus */
  pll_consensus_utree_t * constree = pllmod_utree_weight_consensus(trees,
                                                         weights,
                                                         threshold,
                                                         tree_count);
  if (!constree)
    fatal("ERROR building consensus\n");

  /* print consensus tree in NEWICK format */
  print_newick(constree->tree);

  /* clean up */
  free(weights);
  for (cur_tree=0; cur_tree<tree_count; ++cur_tree)
    pll_utree_destroy(trees[cur_tree], NULL);
  free(trees);

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
    pll_consensus_data_t * cdata = (pll_consensus_data_t *) node->data;
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
