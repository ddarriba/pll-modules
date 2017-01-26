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
#if(USE_HASH)
#include <search.h>
#endif
#include <time.h>

static void fatal (const char * format, ...);

/* function for printing non-binary trees */
static void print_newick_recurse(pll_utree_t * node);
static void print_newick(pll_utree_t * tree);

int main (int argc, char * argv[])
{
  if (argc != 3)
    fatal (" syntax: %s [trees file] [support]", argv[0]);

  double support = atof(argv[2]);

  pll_utree_t * constree =
    pllmod_utree_consensus(argv[1],
                          support,
                          0);
  if (!constree)
    fatal("Error %d: %s\n", pll_errno, pll_errmsg);

  print_newick(constree);

  pll_utree_destroy(constree, NULL);

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

static void print_newick_recurse(pll_utree_t * node)
{
  pll_utree_t * child;
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

static void print_newick(pll_utree_t * tree)
{
  printf("(");
  print_newick_recurse(tree->back);
  pll_utree_t * child = tree->next;
  while(child != tree)
  {
    printf(",");
    print_newick_recurse(child->back);
    child = child->next;
  }
  printf(");\n");
}
