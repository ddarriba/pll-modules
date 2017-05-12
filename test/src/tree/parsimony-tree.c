/*
 Copyright (C) 2016 Alexey Kozlov, Diego Darriba, Tomas Flouri

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

 Contact: Alexey Kozlov <Alexey.Kozlov@h-its.org>,
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

#define ALPHA 0.5
#define N_STATES 4
#define N_CAT_GAMMA 4

#define N_TAXA_BIG 43
#define N_TAXA_SMALL 5
#define N_SITES 491
#define MAX_TAXA 200

#define RAND_SEED 42

int main (int argc, char * argv[])
{
   unsigned int i, n_taxa;
   char * seq[MAX_TAXA], *header[MAX_TAXA];
   long seq_len, header_len, seqno;
   long n_sites = 0;
   pll_fasta_t * fp;
   unsigned int params_indices[4] = {0, 0, 0, 0};
   unsigned int attributes = get_attributes(argc, argv);

   fp = pll_fasta_open ("testdata/small.fas", pll_map_fasta);
   if (!fp)
   {
     printf (" ERROR opening file (%d): %s\n", pll_errno, pll_errmsg);
     exit (PLL_FAILURE);
   }

   /* first read for getting number of taxa and headers */
   i = 0;
   while (pll_fasta_getnext (fp, &header[i], &header_len, &seq[i], &seq_len, &seqno))
   {
     if (!n_sites)
       n_sites = seq_len;
     else if (seq_len != n_sites)
     {
       printf (
           " ERROR: Mismatching sequence length for sequence %d (%ld, and it should be %d)\n",
           i, seq_len - 1, N_SITES);
       exit (PLL_FAILURE);
     }

     printf ("Header of sequence %d(%ld) %s (%ld sites)\n", i, seqno, header[i],
             seq_len);
     ++i;
   }
   pll_fasta_close (fp);
   n_taxa = i;

   if (pll_errno != PLL_ERROR_FILE_EOF)
   {
     printf (" ERROR at the end (%d): %s\n", pll_errno, pll_errmsg);
     exit (PLL_FAILURE);
   }

   unsigned int score;

   pll_utree_t * pars_tree = pllmod_utree_create_parsimony(n_taxa,
                                                      n_sites,
                                                      header,
                                                      seq,
                                                      NULL,  /* site weights */
                                                      pll_map_nt,
                                                      N_STATES,
                                                      attributes,
                                                      RAND_SEED,   /* seed */
                                                      &score);

   pll_unode_t * tree = pars_tree->nodes[2*n_taxa - 3];

   if(!tree)
    fatal("Error creating parsimony [%d]: %s\n", pll_errno, pll_errmsg);

   printf("Parsimony score: %u\n", score);

   pll_utree_show_ascii(tree, PLL_UTREE_SHOW_CLV_INDEX | PLL_UTREE_SHOW_LABEL | PLL_UTREE_SHOW_PMATRIX_INDEX);

    unsigned int tip_nodes_count, inner_nodes_count, nodes_count, branch_count;
    unsigned int matrix_count, ops_count;
    unsigned int * matrix_indices;
    double * branch_lengths;
    pll_partition_t * partition;
    pll_operation_t * operations;
    pll_unode_t ** travbuffer;

    tip_nodes_count = n_taxa;
    inner_nodes_count = n_taxa - 2;
    nodes_count = tip_nodes_count + inner_nodes_count;
    branch_count = 2 * tip_nodes_count - 3;


    printf("SCALE BUFFERS = %d\n", inner_nodes_count);
  /* create the PLL partition instance */
    partition = pll_partition_create (tip_nodes_count,
                                      inner_nodes_count,
                                      N_STATES,
                                      (unsigned int) n_sites,
                                      1,
                                      branch_count,
                                      N_CAT_GAMMA,
                                      inner_nodes_count,
                                      attributes
                                      );

    /* assign tip sequences and free memory */
    for (i=0; i<n_taxa; ++i)
    {
      pll_set_tip_states (partition, i, pll_map_nt, seq[i]);
      free(header[i]);
      free(seq[i]);
    }

    /* initialize base frequencies */
    double frequencies[4] = {0.25,0.25,0.25,0.25};
    pll_set_frequencies (partition, 0, frequencies);

    /* initialize substitution rates */
    double subst_params[6] = {1,1,1,1,1,1};
    pll_set_subst_params (partition, 0,  subst_params);

    /* compute the discretized category rates from a gamma distribution
       with alpha shape 1 and store them in rate_cats  */
    double rate_cats[N_CAT_GAMMA] = { 0 };
    pll_compute_gamma_cats (1, N_CAT_GAMMA, rate_cats, PLL_GAMMA_RATES_MEAN);
    pll_set_category_rates (partition, rate_cats);

    travbuffer = (pll_unode_t **) malloc (nodes_count * sizeof(pll_unode_t *));

    branch_lengths = (double *) malloc (branch_count * sizeof(double));
    matrix_indices = (unsigned int *) malloc (
        branch_count * sizeof(unsigned int));
    operations = (pll_operation_t *) malloc (
        inner_nodes_count * sizeof(pll_operation_t));

    /* perform a postorder traversal of the unrooted tree */
    unsigned int traversal_size;
    if (!pll_utree_traverse (tree,
                             PLL_TREE_TRAVERSE_POSTORDER,
                             cb_full_traversal,
                             travbuffer,
                             &traversal_size))
      fatal ("Function pll_utree_traverse() requires inner nodes as parameters");

    printf("\nTRAVBUFFER: ");
       for (i=0; i<nodes_count; ++i)
         printf("%d/%d/%d  XX ", travbuffer[i]->clv_index, travbuffer[i]->scaler_index, travbuffer[i]->pmatrix_index);
       printf("\n");

    pll_utree_create_operations (travbuffer, traversal_size, branch_lengths,
                                   matrix_indices, operations, &matrix_count,
                                   &ops_count);

    printf("\nMATRICES: ");
    for (i=0; i<matrix_count; ++i)
      printf("%d XX ", matrix_indices[i]);
    printf("\n");
    printf("\nSCALERS: ");
    for (i=0; i<ops_count; ++i)
      printf("%d/%d XX %d/%d XX %d/%d\n", operations[i].child1_clv_index, operations[i].child1_scaler_index, operations[i].child2_clv_index,
             operations[i].child2_scaler_index, operations[i].parent_clv_index, operations[i].parent_scaler_index);
    printf("\n");

    printf ("Traversal size: %d\n", traversal_size);
    printf ("Operations: %d\n", ops_count);
    printf ("Probability Matrices: %d\n", matrix_count);

    pll_update_prob_matrices (partition,
                              params_indices,
                              matrix_indices,
                              branch_lengths,
                              matrix_count);

    pll_update_partials (partition, operations, ops_count);

    double logl = pll_compute_edge_loglikelihood (partition,
                                                  tree->clv_index,
                                                  tree->scaler_index,
                                                  tree->back->clv_index,
                                                  tree->back->scaler_index,
                                                  tree->pmatrix_index,
                                                  params_indices,
                                                  NULL);

   printf ("Log-L at %s-%s: %f\n", tree->label, tree->back->label, logl);

   free (travbuffer);
   free (branch_lengths);
   free (matrix_indices);
   free (operations);

   pll_partition_destroy (partition);
   pll_utree_destroy(pars_tree, NULL);
   return PLL_SUCCESS;
}
