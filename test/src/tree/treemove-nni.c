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
#include "../common.h"

#include <stdarg.h>
#include <search.h>
#include <time.h>

#define STATES    4
#define RATE_CATS 4

#define SHOW_ASCII_TREE 1

#define FASTAFILE "testdata/medium.fas"
#define TREEFILE  "testdata/medium.tree"

typedef struct
{
  int clv_valid;
} node_info_t;

static double evaluate_likelihood(pll_partition_t *partition,
                                   pll_unode_t * tree,
                                   pll_unode_t ** travbuffer,
                                   unsigned int * matrix_indices,
                                   double * branch_lengths,
                                   pll_operation_t * operations,
                                   int update_pmatrices)
{
  double lk;
  unsigned int traversal_size,
               ops_count,
               matrix_count;
  unsigned int params_indices[RATE_CATS] = {0,0,0,0};

  if (update_pmatrices)
  {
    if (!pll_utree_traverse (tree,
                             PLL_TREE_TRAVERSE_POSTORDER,
                             cb_full_traversal,
                             travbuffer,
                             &traversal_size))
      return -1;
  }
  else
  {
    traversal_size = 2;
    travbuffer[0] = tree;
    travbuffer[1] = tree->back;
  }

  pll_utree_create_operations (travbuffer, traversal_size, branch_lengths,
                               matrix_indices, operations, &matrix_count,
                               &ops_count);

  if (update_pmatrices)
  {
    pll_update_prob_matrices (partition, params_indices, matrix_indices, branch_lengths,
                              matrix_count);
  }

  pll_update_partials (partition, operations, ops_count);

  lk = pll_compute_edge_loglikelihood (partition,
                                       tree->clv_index,
                                       tree->scaler_index,
                                       tree->back->clv_index,
                                       tree->back->scaler_index,
                                       tree->pmatrix_index,
                                       params_indices,
                                       NULL);
  return lk;
}

static void apply_move(pll_utree_t * tree,
                       pll_unode_t * edge,
                       int type)
{
  pll_utree_nni(edge, type, 0);

  show_tree(edge, SHOW_ASCII_TREE);

  /* validate tree integrity */
  printf ("Integrity check %s... ", edge->label);
  fflush(stdout);
  pll_errno = 0;
  if (!pll_utree_check_integrity (tree))
    fatal ("Tree is not consistent %ld %s", pll_errno, pll_errmsg);
  printf ("OK!\n");
}

int main (int argc, char *argv[])
{
  double logl, logl_start, logl_end;

  unsigned int i;
  unsigned int tip_nodes_count, inner_nodes_count, nodes_count, branch_count;
  unsigned int * matrix_indices;
  double * branch_lengths;
  pll_partition_t * partition;
  pll_operation_t * operations;
  pll_unode_t ** travbuffer;

  unsigned int attributes = get_attributes(argc, argv);

  /* parse the unrooted binary tree in newick format, and store the number
   of tip nodes in tip_nodes_count */
  printf("Parsing tree: %s\n", TREEFILE);
  pll_utree_t * parsed_tree = pll_utree_parse_newick (TREEFILE);
  pll_unode_t * tree = 0;
  if (!parsed_tree)
    fatal ("Error parsing %s [%d] %s", TREEFILE, pll_errno, pll_errmsg);
  tip_nodes_count = parsed_tree->tip_count;

  /* compute and show node count information */
  inner_nodes_count = tip_nodes_count - 2;
  nodes_count = inner_nodes_count + tip_nodes_count;
  branch_count = nodes_count - 1;

  printf ("  Number of tip/leaf nodes in tree: %d\n", tip_nodes_count);
  printf ("  Number of inner nodes in tree: %d\n", inner_nodes_count);
  printf ("  Total number of nodes in tree: %d\n", nodes_count);
  printf ("  Number of branches in tree: %d\n", branch_count);

  /*  obtain an array of pointers to tip and inner nodes */
  pll_unode_t ** tipnodes = parsed_tree->nodes;
  pll_unode_t ** innernodes = parsed_tree->nodes + tip_nodes_count;

  /* place the virtual root at a random inner node */
  tree = innernodes[(unsigned int) rand() % inner_nodes_count];

  show_tree (tree, SHOW_ASCII_TREE);

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
  printf("Reading FASTA file: %s\n", FASTAFILE);
  pll_fasta_t * fp = pll_fasta_open (FASTAFILE, pll_map_fasta);
  if (!fp)
    fatal ("%s does not exist", FASTAFILE);

  char * seq = NULL;
  char * hdr = NULL;
  long seqlen;
  long hdrlen;
  long seqno;

  /* allocate arrays to store FASTA headers and sequences */
  char ** headers = (char **) calloc (tip_nodes_count, sizeof(char *));
  char ** seqdata = (char **) calloc (tip_nodes_count, sizeof(char *));

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
  double frequencies[4] = {0.25,0.25,0.25,0.25};
  pll_set_frequencies (partition, 0, frequencies);

  /* initialize substitution rates */
  double subst_params[6] = {1,1,1,1,1,1};
  pll_set_subst_params (partition, 0, subst_params);

  /* compute the discretized category rates from a gamma distribution
     with alpha shape 1 and store them in rate_cats  */
  double rate_cats[RATE_CATS] = { 0 };
  pll_compute_gamma_cats (1, RATE_CATS, rate_cats, PLL_GAMMA_RATES_MEAN);
  pll_set_category_rates (partition, rate_cats);

  printf ("Model paramters:\n");
  printf ("  Frequencies:  ");
  for (i=0; i<STATES; i++)
    printf("%.4f ", partition->frequencies[0][i]);
  printf("\n");
  printf ("  Subst. rates: ");
  for (i=0; i<(STATES*(STATES-1))/2; i++)
    printf("%.4f ", partition->subst_params[0][i]);
  printf("\n");
  printf ("  Gamma rates:  ");
  for (i=0; i<RATE_CATS; i++)
    printf("%.4f ", partition->rates[i]);
  printf("\n");

  /* allocate a buffer for storing pointers to nodes of the tree in postorder
   traversal */
  travbuffer = (pll_unode_t **) malloc (nodes_count * sizeof(pll_unode_t *));

  branch_lengths = (double *) malloc (branch_count * sizeof(double));
  matrix_indices = (unsigned int *) malloc (
      branch_count * sizeof(unsigned int));
  operations = (pll_operation_t *) malloc (
      inner_nodes_count * sizeof(pll_operation_t));

  evaluate_likelihood(partition, tree, travbuffer,
                      matrix_indices, branch_lengths,
                      operations, 1);

  for (i=0; i<5; i++)
  {
    printf("Iteration %d\n", i);

    logl_start = logl = evaluate_likelihood(partition, tree, travbuffer,
                        matrix_indices, branch_lengths,
                        operations, 1);
    printf ("Log-L[ST] at %s-%s: %f\n", tree->label, tree->back->label, logl);

    /* Test NNI */
    pll_unode_t * nni_edge;

    printf("\n\n");
    printf("Obtaining random inner edge\n");
    nni_edge = innernodes[(unsigned int) rand () % inner_nodes_count];

    show_tree(nni_edge, SHOW_ASCII_TREE);

    /* update CLVs towards new root */
    logl = evaluate_likelihood(partition, nni_edge, travbuffer,
                        matrix_indices, branch_lengths,
                        operations, 1);

    assert(fabs(logl - logl_start) < 1e-7);

    printf ("Tree move at %s-%s\n", nni_edge->label, nni_edge->back->label);

    apply_move(parsed_tree, nni_edge, PLL_UTREE_MOVE_NNI_LEFT);

    logl = evaluate_likelihood(partition, nni_edge, travbuffer,
                        matrix_indices, branch_lengths,
                        operations, 0);
    printf ("Log-L[M1] at %s-%s: %f\n", nni_edge->label, nni_edge->back->label, logl);

    /* second move */
    apply_move(parsed_tree, nni_edge, PLL_UTREE_MOVE_NNI_RIGHT);

    logl = evaluate_likelihood(partition, nni_edge, travbuffer,
                        matrix_indices, branch_lengths,
                        operations, 0);
    printf ("Log-L[M2] at %s-%s: %f\n", nni_edge->label, nni_edge->back->label, logl);

    /* rollback */
    apply_move(parsed_tree, nni_edge, PLL_UTREE_MOVE_NNI_LEFT);

    logl_end = logl = evaluate_likelihood(partition, nni_edge, travbuffer,
                        matrix_indices, branch_lengths,
                        operations, 0);
    printf ("Log-L[RB] at %s-%s: %f\n", nni_edge->label, nni_edge->back->label, logl);

    if (fabs(logl_start - logl_end) > 1e-7)
    {
      printf("Error: Starting and final Log-LK differ\n");
      assert(0);
    }
  }

  /* clean */
  pll_partition_destroy (partition);
  free (travbuffer);
  free (branch_lengths);
  free (matrix_indices);
  free (operations);
  pll_utree_destroy (parsed_tree, NULL);

  printf("Test OK!\n");

  return (EXIT_SUCCESS);
}
