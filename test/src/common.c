#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

unsigned int get_attributes(int argc, char **argv)
{
  int i;
  unsigned int attributes = PLL_ATTRIB_ARCH_CPU;

  for (i=1; i<argc; ++i)
  {
    if (!strcmp (argv[i], "tv"))
    {
      /* tipvector */
      attributes |= PLL_ATTRIB_PATTERN_TIP;
    }
    else if (!strcmp (argv[i], "avx"))
    {
      /* avx vectorization */
      attributes |= PLL_ATTRIB_ARCH_AVX;
    }
    else if (!strcmp (argv[i], "sse"))
    {
      /* sse3 vectorization */
      attributes |= PLL_ATTRIB_ARCH_SSE;
    }
    else if (!strcmp (argv[i], "sr"))
    {
      attributes |= PLL_ATTRIB_SITE_REPEATS;
    }
    else
    {
      printf("Unrecognised attribute: %s\n", argv[i]);
      exit(1);
    }
  }
    return attributes;
}

void skip_test ()
{
  printf ("Skip\n");
  exit (0);
}


int cb_full_traversal (pll_unode_t * node)
{
  return 1;
}

int cb_rfull_traversal (pll_rnode_t * node)
{
  return 1;
}

void show_tree (pll_unode_t * tree, int SHOW_ASCII_TREE)
{
  if(SHOW_ASCII_TREE)
  {
    printf ("\n");
    pll_utree_show_ascii (
        tree,
        PLL_UTREE_SHOW_LABEL |
        PLL_UTREE_SHOW_BRANCH_LENGTH |
        PLL_UTREE_SHOW_CLV_INDEX | PLL_UTREE_SHOW_PMATRIX_INDEX
            | PLL_UTREE_SHOW_SCALER_INDEX);
    char * newick = pll_utree_export_newick (tree, NULL);
    printf ("%s\n\n", newick);
    free (newick);
  }
  else
  {
    printf ("ASCII tree not shown (SHOW_ASCII_TREE flag)\n");
    return;
  }
}

void show_rtree (pll_rnode_t * tree, int SHOW_ASCII_TREE)
{
  if(SHOW_ASCII_TREE)
  {
    printf ("\n");
    pll_rtree_show_ascii (
        tree,
        PLL_UTREE_SHOW_LABEL |
        PLL_UTREE_SHOW_BRANCH_LENGTH |
        PLL_UTREE_SHOW_CLV_INDEX | PLL_UTREE_SHOW_PMATRIX_INDEX
            | PLL_UTREE_SHOW_SCALER_INDEX);
    char * newick = pll_rtree_export_newick (tree, NULL);
    printf ("%s\n\n", newick);
    free (newick);
  }
  else
  {
    printf ("ASCII tree not shown (SHOW_ASCII_TREE flag)\n");
    return;
  }
}

void fatal(const char * format, ...)
{
  va_list argptr;
  va_start(argptr, format);
  vfprintf(stderr, format, argptr);
  va_end(argptr);
  fprintf(stderr, "\n");
  exit(EXIT_FAILURE);
}
