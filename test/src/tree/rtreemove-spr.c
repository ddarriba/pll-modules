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
#include "../common.h"

#include <stdarg.h>
#include <search.h>
#include <time.h>

#define STATES    4
#define RATE_CATS 4

#define SHOW_ASCII_TREE 1

#define FASTAFILE "testdata/medium.fas"
#define TREEFILE  "testdata/medium.rooted.tree"

typedef struct
{
  int clv_valid;
} node_info_t;

static double evaluate_likelihood (pll_partition_t *partition,
                                   pll_rnode_t * tree,
                                   pll_rnode_t ** travbuffer,
                                   unsigned int * matrix_indices,
                                   double * branch_lengths,
                                   pll_operation_t * operations)
{
  double lk;
  unsigned int traversal_size;
  unsigned int ops_count;
  unsigned int matrix_count;
  unsigned int params_indices[RATE_CATS];
  int i;

  for (i=0; i<RATE_CATS; ++i)
    params_indices[i] = 0;

  if (!pll_rtree_traverse (tree, PLL_TREE_TRAVERSE_POSTORDER,
                           cb_rfull_traversal, travbuffer,
                           &traversal_size))
    return -1;

  pll_rtree_create_operations (travbuffer, traversal_size, branch_lengths,
                               matrix_indices, operations, &matrix_count,
                               &ops_count);

  pll_update_prob_matrices (partition, params_indices, matrix_indices,
                            branch_lengths, matrix_count);

  pll_update_partials (partition, operations, ops_count);

  lk = pll_compute_root_loglikelihood (partition, tree->clv_index,
                                       tree->scaler_index,
                                       params_indices,
                                       NULL);
  return lk;
}

static void apply_move (pll_rnode_t * edge, pll_rnode_t * tree,
                        pll_rnode_t ** root,
                        pll_tree_rollback_t * rollback_stack,
                        int * rollback_stack_top)
{
  printf ("Apply SPR: %s to %s, with root %s\n", edge->label, tree->label, (*root)->label);
  if (!pllmod_rtree_spr (edge, tree, root, &rollback_stack[++(*rollback_stack_top)]))
  {
    printf ("Error %d: %s\n", pll_errno, pll_errmsg);
    exit (1);
  }

  show_rtree (*root, SHOW_ASCII_TREE);
}

static void undo_move (pll_rnode_t **root,
                       pll_tree_rollback_t * rollback_stack,
                       int * rollback_stack_top)
{
  assert(!(*root)->parent);

  printf("Rollback move!\n");
  pllmod_tree_rollback (&rollback_stack[(*rollback_stack_top)--]);

  while((*root)->parent) *root = (*root)->parent;
}

#define ROLLBACK_STACK_SIZE 10
#define EMPTY_STACK         -1

int main (int argc, char * argv[])
{
  double logl, logl_start, logl_end;

  unsigned int i;
  unsigned int tip_nodes_count, inner_nodes_count, nodes_count, branch_count;
  unsigned int * matrix_indices;
  double * branch_lengths;
  pll_partition_t * partition;
  pll_operation_t * operations;
  pll_rnode_t ** travbuffer;

  char * seq = NULL;
  char * hdr = NULL;
  long seqlen;
  long hdrlen;
  long seqno;
  char ** headers;
  char ** seqdata;
  pll_tree_rollback_t * rollback_stack = (pll_tree_rollback_t *) malloc (
      ROLLBACK_STACK_SIZE * sizeof(pll_tree_rollback_t));
  int rollback_stack_top = EMPTY_STACK;

  unsigned int attributes = get_attributes(argc, argv);

  /* parse the unrooted binary tree in newick format, and store the number
   of tip nodes in tip_nodes_count */
  printf ("Parsing tree: %s\n", TREEFILE);
  pll_rtree_t * parsed_tree = pll_rtree_parse_newick (TREEFILE);
  if (!parsed_tree)
    fatal ("Error parsing %s", TREEFILE);
  tip_nodes_count = parsed_tree->tip_count;
  pll_rnode_t * tree = parsed_tree->root;

  /* compute and show node count information */
  inner_nodes_count = tip_nodes_count - 1;
  nodes_count = inner_nodes_count + tip_nodes_count;
  branch_count = nodes_count - 1;

  printf ("  Number of tip/leaf nodes in tree: %d\n", tip_nodes_count);
  printf ("  Number of inner nodes in tree: %d\n", inner_nodes_count);
  printf ("  Total number of nodes in tree: %d\n", nodes_count);
  printf ("  Number of branches in tree: %d\n", branch_count);

  /*  obtain an array of pointers to tip and inner nodes */
  pll_rnode_t ** tipnodes = parsed_tree->nodes;
  pll_rnode_t ** innernodes = parsed_tree->nodes + tip_nodes_count;

  show_rtree (tree, SHOW_ASCII_TREE);

  /* create a libc hash table of size tip_nodes_count */
  hcreate (tip_nodes_count);

  /* populate a libc hash table with tree tip labels */
  unsigned int * data = (unsigned int *) malloc (
      tip_nodes_count * sizeof(unsigned int));
  for (i = 0; i < tip_nodes_count; ++i)
  {
    data[i] = i;
    ENTRY entry;
    entry.key = tipnodes[i]->label;
    entry.data = (void *) (data + i);
    hsearch (entry, ENTER);
  }

  /* open FASTA file */
  printf ("Reading FASTA file: %s\n", FASTAFILE);
  pll_fasta_t * fp = pll_fasta_open (FASTAFILE, pll_map_fasta);
  if (!fp)
    fatal ("%s does not exist", FASTAFILE);

  /* allocate arrays to store FASTA headers and sequences */
  headers = (char **) calloc (tip_nodes_count, sizeof(char *));
  seqdata = (char **) calloc (tip_nodes_count, sizeof(char *));

  /* read FASTA sequences and make sure they are all of the same length */
  int sites = -1;
  for (i = 0; pll_fasta_getnext (fp, &hdr, &hdrlen, &seq, &seqlen, &seqno); ++i)
  {
    if (i >= tip_nodes_count)
      fatal ("FASTA file contains more sequences than expected");

    if (sites != -1 && sites != seqlen)
      fatal ("FASTA file does not contain equal size sequences\n");

    if (sites == -1)
      sites = (int) seqlen;

    headers[i] = hdr;
    seqdata[i] = seq;
  }

  /* did we stop reading the file because we reached EOF? */
  if (pll_errno != PLL_ERROR_FILE_EOF)
    fatal ("Error in %s", FASTAFILE);

  /* close FASTA file */
  pll_fasta_close (fp);

  if (sites == -1)
    fatal ("Unable to read alignment");

  if (i != tip_nodes_count)
    fatal ("Some taxa are missing from FASTA file");

  printf ("  Length of sequences: %d\n", sites);

  /* create the PLL partition instance */
  partition = pll_partition_create (tip_nodes_count,
                                    inner_nodes_count,
                                    STATES,
                                    (unsigned int) sites,
                                    1,
                                    branch_count,
                                    RATE_CATS,
                                    inner_nodes_count,
                                    attributes
                                    );

  /* find sequences in hash table and link them with the corresponding taxa */
  for (i = 0; i < tip_nodes_count; ++i)
  {
    ENTRY query;
    query.key = headers[i];
    ENTRY * found = NULL;

    found = hsearch (query, FIND);

    if (!found)
      fatal ("Sequence with header %s does not appear in the tree", hdr);

    unsigned int tip_clv_index = *((unsigned int *) (found->data));

    pll_set_tip_states (partition, tip_clv_index, pll_map_nt, seqdata[i]);
  }

  /* destroy hash table */
  hdestroy ();

  /* we no longer need these two arrays (keys and values of hash table... */
  free (data);

  /* ...neither the sequences and the headers as they are already
   present in the form of probabilities in the tip CLVs */
  for (i = 0; i < tip_nodes_count; ++i)
  {
    free (seqdata[i]);
    free (headers[i]);
  }
  free (seqdata);
  free (headers);

  /* initialize base frequencies */
  double frequencies[4] =
    { 0.25, 0.25, 0.25, 0.25 };
  pll_set_frequencies (partition, 0, frequencies);

  /* initialize substitution rates */
  double subst_params[6] =
    { 1, 1, 1, 1, 1, 1 };
  pll_set_subst_params (partition, 0, subst_params);

  /* compute the discretized category rates from a gamma distribution
   with alpha shape 1 and store them in rate_cats  */
  double rate_cats[RATE_CATS] =
    { 0 };
  pll_compute_gamma_cats (1, RATE_CATS, rate_cats, PLL_GAMMA_RATES_MEAN);
  pll_set_category_rates (partition, rate_cats);

  printf ("Model paramters:\n");
  printf ("  Frequencies:  ");
  for (i = 0; i < STATES; i++)
    printf ("%.4f ", partition->frequencies[0][i]);
  printf ("\n");
  printf ("  Subst. rates: ");
  for (i = 0; i < (STATES * (STATES - 1)) / 2; i++)
    printf ("%.4f ", partition->subst_params[0][i]);
  printf ("\n");
  printf ("  Gamma rates:  ");
  for (i = 0; i < RATE_CATS; i++)
    printf ("%.4f ", partition->rates[i]);
  printf ("\n");

  /* allocate a buffer for storing pointers to nodes of the tree in postorder
   traversal */
  travbuffer = (pll_rnode_t **) malloc (nodes_count * sizeof(pll_rnode_t *));

  branch_lengths = (double *) malloc (branch_count * sizeof(double));
  matrix_indices = (unsigned int *) malloc (
      branch_count * sizeof(unsigned int));
  operations = (pll_operation_t *) malloc (
      inner_nodes_count * sizeof(pll_operation_t));

  unsigned int distance = 8;
  unsigned int n_nodes_at_dist;
  pll_rnode_t ** nodes_at_dist = (pll_rnode_t **) malloc (
      sizeof(pll_rnode_t *) * pow (2, distance + 1));

  for (i = 0; i < 5; i++)
  {
    printf ("Iteration %d\n", i);

    logl_start = logl = evaluate_likelihood (partition, tree, travbuffer,
                                             matrix_indices, branch_lengths,
                                             operations);
    show_rtree (tree, SHOW_ASCII_TREE);
    printf ("Log-L[ST] %f\n", logl);

    /* Test SPR */
    pll_rnode_t *prune_edge, *regraft_edge;

    printf ("\n\n");
    printf ("Obtaining random inner edge\n");

    do {
      prune_edge = innernodes[(unsigned int) rand () % inner_nodes_count];
    } while(!prune_edge->parent);
    printf ("Tree prune at %s\n", prune_edge->label);

    pllmod_rtree_nodes_at_node_dist (prune_edge, nodes_at_dist, &n_nodes_at_dist,
                                     2, distance);
    assert(n_nodes_at_dist > 0);

    int j;
    for (j=0;j<n_nodes_at_dist;++j)
      printf("%s ", nodes_at_dist[j]->label);
    printf("\n");

    regraft_edge = nodes_at_dist[(unsigned int) rand () % n_nodes_at_dist];

    apply_move (prune_edge, regraft_edge, &tree, rollback_stack, &rollback_stack_top);

    logl = evaluate_likelihood (partition, tree, travbuffer,
                                matrix_indices, branch_lengths, operations);

    printf ("Log-L[M1] at %s: %f\n", tree->label, logl);

    printf ("Obtaining random inner edge\n");
    do
    {
      prune_edge = innernodes[(unsigned int) rand () % inner_nodes_count];
    } while(!prune_edge->parent);

    pllmod_rtree_nodes_at_node_dist (prune_edge, nodes_at_dist, &n_nodes_at_dist,
                                     2, distance);
    assert(n_nodes_at_dist > 0);

    regraft_edge = nodes_at_dist[(unsigned int) rand () % n_nodes_at_dist];

    printf ("Tree prune at %s\n", prune_edge->label);
    printf ("Tree regraft at %s\n", regraft_edge->label);

    apply_move (prune_edge, regraft_edge, &tree, rollback_stack, &rollback_stack_top);

    logl = evaluate_likelihood (partition, tree, travbuffer,
                                matrix_indices, branch_lengths, operations);
    printf ("Log-L[M1] at %s: %f\n", tree->label, logl);

    /* rollback */
    printf ("Tree prune at %s\n", prune_edge->label);
    printf ("Tree regraft at %s\n", regraft_edge->label);

    /* rollback */

    /* undo first move */
    undo_move (&tree, rollback_stack, &rollback_stack_top);
    /* undo second move */
    undo_move (&tree, rollback_stack, &rollback_stack_top);

    show_rtree (prune_edge, SHOW_ASCII_TREE);

    logl_end = logl = evaluate_likelihood (partition, tree, travbuffer,
                                           matrix_indices, branch_lengths,
                                           operations);
    printf ("Log-L[RB] at %s: %f\n", tree->label, logl);

    if (fabs (logl_start - logl_end) > 1e-7)
    {
      printf ("Error: Starting and final Log-LK differ\n");
      assert(0);
    }
  }

  assert (rollback_stack_top == EMPTY_STACK);

  /* clean */
  pll_partition_destroy (partition);
  free (rollback_stack);
  free (nodes_at_dist);
  free (travbuffer);
  free (branch_lengths);
  free (matrix_indices);
  free (operations);
  pll_rtree_destroy (parsed_tree, NULL);

  printf ("Test OK!\n");

  return (EXIT_SUCCESS);
}
