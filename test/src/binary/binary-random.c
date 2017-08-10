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
#include "pll_binary.h"
#include "../common.h"

#include <string.h>
#include <search.h>

#define ALPHA        0.5
#define N_STATES      4
#define N_SUBST_RATES 6
#define N_RATE_CATS   4

#define MSA_FILENAME  "testdata/246x4465.fas"
#define TREE_FILENAME "testdata/246x4465.tree"
// #define MSA_FILENAME  "testdata/small.fas"
// #define TREE_FILENAME "testdata/small.tree"

#define BLOCK_ID_PARTITION 1000
#define BLOCK_ID_TREE      2000
#define BLOCK_ID_CLV       3000
#define BLOCK_ID_CUSTOM    4000
#define BLOCK_ID_REPEATS   5000

static int cb_traversal(pll_unode_t * node)
{
  return 1;
}

static pll_partition_t * parse_msa(unsigned int attributes, pll_utree_t * tree)
{
  unsigned int i;
  unsigned int tip_clv_index;
  char * seq, *header;
  char ** headers, ** seqdata;
  long seq_len    = 0,
       read_len   = 0,
       header_len = 0,
       seqno      = 0;
  pll_fasta_t * fp;
  pll_partition_t * partition;

  pll_unode_t ** tipnodes = tree->nodes;
  unsigned int tip_nodes_count = tree->tip_count;

  headers = (char **)calloc(tip_nodes_count, sizeof(char *));
  seqdata = (char **)calloc(tip_nodes_count, sizeof(char *));

  /* create a libc hash table of size tip_nodes_count */
  hcreate(tip_nodes_count);

  /* populate a libc hash table with tree tip labels */
  unsigned int * data = (unsigned int *) malloc ( tip_nodes_count *
                                                  sizeof(unsigned int) );
  for (i = 0; i < tip_nodes_count; ++i)
  {
    data[i] = i;
    ENTRY entry;
    entry.key = tipnodes[i]->label;
    entry.data = (void *)(data+i);
    hsearch(entry, ENTER);
  }

  fp = pll_fasta_open (MSA_FILENAME, pll_map_fasta);
  if (!fp)
  {
    printf (" ERROR opening file (%d): %s\n", pll_errno, pll_errmsg);
    return NULL;
  }

  i = 0;
  while (pll_fasta_getnext (fp, &header, &header_len, &seq, &read_len, &seqno))
  {
    if (!seq_len)
    {
      seq_len = read_len;
    }
    else if (seq_len != read_len)
    {
      printf (
          " ERROR: Mismatching sequence length for sequence %d (%ld, and it should be %ld)\n",
          i, read_len - 1, seq_len);
      return NULL;
    }
    headers[i] = header;
    seqdata[i] = seq;
    ++i;
  }

  if (pll_errno != PLL_ERROR_FILE_EOF)
  {
    printf (" ERROR at the end (%d): %s\n", pll_errno, pll_errmsg);
    return NULL;
  }

  pll_fasta_close (fp);

  partition = pll_partition_create(tip_nodes_count,      /* tips */
                                    tip_nodes_count - 2, /* clv buffers */
                                    N_STATES,            /* states */
                                    seq_len,             /* sites */
                                    1,           /* different rate parameters */
                                    2*tip_nodes_count - 3, /* prob matrices */
                                    N_RATE_CATS,           /* rate categories */
                                    tip_nodes_count - 2,   /* scale buffers */
                                    attributes             /* attributes */
                                    );

  for (i = 0; i < tip_nodes_count; ++i)
  {
    ENTRY query;
    query.key = headers[i];
    ENTRY * found = NULL;

    found = hsearch(query,FIND);

    if (!found)
      fatal("Sequence with header %s does not appear in the tree", headers[i]);

    tip_clv_index = *((unsigned int *)(found->data));
    pll_set_tip_states(partition, tip_clv_index, pll_map_nt, seqdata[i]);
  }

  hdestroy();
  free(data);

  for(i = 0; i < tip_nodes_count; ++i)
  {
    free(seqdata[i]);
    free(headers[i]);
  }
  free(seqdata);
  free(headers);

  return partition;
}

typedef struct
{
  int * a;
  int * b;
} test_t;

int write(void * data, size_t size, size_t count, FILE * file)
{
  size_t ret = fwrite(data, size, count, file);
  if (ret != count)
  {
    return PLL_FAILURE;
  }
  return PLL_SUCCESS;
}

int main (int argc, char * argv[])
{
  unsigned int attributes = get_attributes(argc, argv);
  unsigned int matrix_count, ops_count;
  unsigned int tip_nodes_count, inner_nodes_count, nodes_count, branch_count;
  pll_partition_t * partition;
  pll_unode_t * tree;
  double logl, save_logl;
  int i;

  pll_unode_t ** travbuffer;
  double * branch_lengths;
  unsigned int * matrix_indices;
  pll_operation_t * operations;

  double frequencies[N_STATES]             = { 0.17, 0.19, 0.25, 0.39 };
  double subst_params[N_SUBST_RATES]       = {1,1,1,1,1,1};
  double rate_cats[N_RATE_CATS]            = {0};
  unsigned int params_indices[N_RATE_CATS] = {0, 0, 0, 0};

  unsigned int lk_parent_clv_index, lk_parent_scaler_index,
               lk_child_clv_index, lk_child_scaler_index,
               lk_pmatrix_index;

  pll_utree_t * parsed_tree = pll_utree_parse_newick(TREE_FILENAME);
  tip_nodes_count = parsed_tree->tip_count;
  tree = parsed_tree->nodes[2*tip_nodes_count - 3];

  if (!tree)
  {
    printf ("Error %d parsing tree: %s\n", pll_errno, pll_errmsg);
    return 1;
  }

  inner_nodes_count = tip_nodes_count - 2;
  nodes_count = inner_nodes_count + tip_nodes_count;
  branch_count = nodes_count - 1;

  printf("Number of tip/leaf nodes in tree: %d\n", tip_nodes_count);
  printf("Number of inner nodes in tree: %d\n", inner_nodes_count);
  printf("Total number of nodes in tree: %d\n", nodes_count);
  printf("Number of branches in tree: %d\n", branch_count);

  partition = parse_msa (attributes, parsed_tree);
  if (!partition)
  {
    printf ("Error creating partition\n");
    return 1;
  }

  pll_compute_gamma_cats(ALPHA, N_RATE_CATS, rate_cats, PLL_GAMMA_RATES_MEAN);
  pll_set_frequencies(partition, 0, frequencies);
  pll_set_subst_params(partition, 0, subst_params);
  pll_set_category_rates(partition, rate_cats);

  travbuffer = (pll_unode_t **)malloc(nodes_count * sizeof(pll_unode_t *));
  branch_lengths = (double *)malloc(branch_count * sizeof(double));
  matrix_indices = (unsigned int *)malloc(branch_count * sizeof(int));
  operations = (pll_operation_t *)malloc(inner_nodes_count *
                                                sizeof(pll_operation_t));

  unsigned int traversal_size;

  pll_unode_t * node = tree;

  if (!pll_utree_traverse(node,
                          PLL_TREE_TRAVERSE_POSTORDER,
                          cb_traversal,
                          travbuffer,
                          &traversal_size))
    fatal("Function pll_utree_traverse() requires inner nodes as parameters");

  pll_utree_create_operations(travbuffer,
                              traversal_size,
                              branch_lengths,
                              matrix_indices,
                              operations,
                              &matrix_count,
                              &ops_count);

  printf("\nComputing logL between CLV %d and %d - "
         "(pmatrix %d with branch length %f)\n",
            node->clv_index,
            node->back->clv_index,
            node->pmatrix_index,
            node->length);

  printf ("Traversal size: %d\n", traversal_size);
  printf ("Operations: %d\n", ops_count);
  printf ("Matrices: %d\n", matrix_count);

  pll_update_prob_matrices(partition,
                           params_indices,
                           matrix_indices,
                           branch_lengths,
                           matrix_count);

  pll_update_partials(partition, operations, ops_count);

  lk_parent_clv_index = node->clv_index;
  lk_parent_scaler_index = node->scaler_index;
  lk_child_clv_index = node->back->clv_index;
  lk_child_scaler_index = node->back->scaler_index;
  lk_pmatrix_index = node->pmatrix_index;

  logl = pll_compute_edge_loglikelihood(partition,
                                        lk_parent_clv_index,
                                        lk_parent_scaler_index,
                                        lk_child_clv_index,
                                        lk_child_scaler_index,
                                        lk_pmatrix_index,
                                        params_indices,
                                        NULL);

  save_logl = logl;
  printf("Log-L: %f\n", logl);

  FILE * bin_file;
  pll_binary_header_t bin_header;
  const char * bin_fname = "test.bin";

  printf("** create binary file\n");
  bin_file = pllmod_binary_create(bin_fname,
                               &bin_header,
                               PLLMOD_BIN_ACCESS_RANDOM,
                               10); /* allocate for up to 10 blocks */

  if (!bin_file)
  {
    fatal("Cannot create binary file: %s\n", bin_fname);
  }

  printf("** dump partition\n");

  /* We save the structures in an arbitrary order */

  /* IMPORTANT! Attribute PLLMOD_BIN_ATTRIB_UPDATE_MAP must be set! */
  if (!pllmod_binary_partition_dump(bin_file,
                            BLOCK_ID_PARTITION,
                            partition,
                            PLLMOD_BIN_ATTRIB_PARTITION_DUMP_CLV |
                            PLLMOD_BIN_ATTRIB_PARTITION_DUMP_WGT |
                            PLLMOD_BIN_ATTRIB_UPDATE_MAP))
  {
    printf("Error dumping partition\n");
  }

  /* dump tree */
  if (!pllmod_binary_utree_dump(bin_file,
                       BLOCK_ID_TREE,
                       tree,
                       tip_nodes_count,
                       PLLMOD_BIN_ATTRIB_UPDATE_MAP))
  {
    printf("Error dumping tree\n");
  }

  /* dump 5 arbitrary CLVs and save original values*/
  const int n_clvs   = 5;
  size_t clv_size = partition->states_padded * partition->rate_cats * partition->sites;

  double * saved_clvs[n_clvs];
  for (unsigned int i = 0; i < n_clvs; ++i)
  {
    saved_clvs[i] = (double *) malloc(pll_get_clv_size(partition, partition->tips +  i) * sizeof(double));
  }
  for (i=0; i<n_clvs; ++i)
  {
    int clv_index;
    clv_index = partition->tips + i;
    assert (clv_index < (partition->tips + partition->clv_buffers));
    pllmod_binary_clv_dump(bin_file,
                        BLOCK_ID_CLV + i,
                        partition,
                        clv_index,
                        PLLMOD_BIN_ATTRIB_UPDATE_MAP);
    memcpy(saved_clvs[i],
           partition->clv[clv_index],
           sizeof(double) * pll_get_clv_size(partition, clv_index));
  }

  printf("** close binary file\n");

  pllmod_binary_close(bin_file);

  // pll_utree_show_ascii(tree, (1<<5)-1);

  /* clean */
  pll_partition_destroy(partition);
  pll_utree_destroy(parsed_tree, NULL);

  printf("\n\n");


  /* reload */
  printf("** reload data\n");
  pll_binary_header_t input_header;
  unsigned int bin_attributes = 0;
  pll_block_map_t * block_map;
  unsigned int n_blocks;

  bin_file = pllmod_binary_open(bin_fname, &input_header);

  block_map = pllmod_binary_get_map(bin_file, &n_blocks);

  printf("There are %d blocks in the map\n", n_blocks);
  int partition_offset = 0;
  for (i=0; i<n_blocks; ++i)
  {
    /* we cannot print the offset for the test since PATTERN_TIP
    affects the size of the CLVs */
    printf("%5ld\n", block_map[i].block_id);
    // printf("%5ld %20ld\n", block_map[i].block_id, block_map[i].block_offset);
    //
    if (block_map[i].block_id == BLOCK_ID_PARTITION)
    {
      partition_offset = block_map[i].block_offset;
    }
  }

  /* For the offset we can use the actual offset (from the block_map),
     or PLLMOD_BIN_ACCESS_SEEK */
  partition = pllmod_binary_partition_load(bin_file,
                                        BLOCK_ID_PARTITION,
                                        NULL, /* in order to create a new partition */
                                        &bin_attributes,
                                        partition_offset);

  if (!partition)
    printf("Error loading partition\n");

  double * loaded_clvs[n_clvs];
  for (unsigned int i = 0; i < n_clvs; ++i)
  {
    loaded_clvs[i] = (double *) malloc(pll_get_clv_size(partition, partition->tips + i) * sizeof(double));
  }
  /* load clvs in reverse order for testing */
  for (i=(n_clvs-1); i>=0; --i)
  {
    int clv_index = partition->tips + i;
    assert (clv_index < (partition->tips + partition->clv_buffers));

    /* reset involved clvs */
    memset(partition->clv[clv_index], 0, sizeof(double) * pll_get_clv_size(partition, clv_index));

    if (!pllmod_binary_clv_load(bin_file,
                        BLOCK_ID_CLV + i,
                        partition,
                        clv_index,
                        &bin_attributes,
                        PLLMOD_BIN_ACCESS_SEEK))
    {
      printf("Error loading CLV %d\n", clv_index);
      printf("%d : %s\n", pll_errno, pll_errmsg);
      exit(1);
    }

    memcpy(loaded_clvs[i],
           partition->clv[clv_index],
           sizeof(double) * pll_get_clv_size(partition, clv_index));
    // loaded_clvs_ptr += clv_size;
  }

  /* compare CLVs */
  for (unsigned int i = 0; i < n_clvs; ++i) 
  {
    if (memcmp(saved_clvs[i], loaded_clvs[i], pll_get_clv_size(partition, partition->tips + i) * sizeof(double)))
    {
      printf("Error! CLVs do not agree\n");
      exit(1);
    }
  }
  printf("saved and restored clvs OK!\n");

  for (unsigned int i = 0; i < n_clvs; ++i) {
    free(saved_clvs[i]);
    free(loaded_clvs[i]);
  }

  /* first check if partition was restored correctly */
  logl = pll_compute_edge_loglikelihood(partition,
                                         lk_parent_clv_index,
                                         lk_parent_scaler_index,
                                         lk_child_clv_index,
                                         lk_child_scaler_index,
                                         lk_pmatrix_index,
                                         params_indices,
                                         NULL);

   printf("Restored Log-L: %f\n", logl);
   if (fabs(logl - save_logl) < 1e-7)
     printf("Likelihoods OK!\n");
   else
     fatal("Error: Saved and loaded logL do not agree!! %f %f\n", logl, save_logl);

  /* new we try with ACCESS_SEEK instead of the value taken from the map */
  tree = pllmod_binary_utree_load(bin_file,
                               BLOCK_ID_TREE,
                               &bin_attributes,
                               PLLMOD_BIN_ACCESS_SEEK);
  if (!tree)
    fatal("Error loading tree!\n");

pllmod_binary_close(bin_file);

  if (!pll_utree_traverse(tree,
                          PLL_TREE_TRAVERSE_POSTORDER,
                          cb_traversal,
                          travbuffer,
                          &traversal_size))
    fatal("Function pll_utree_traverse() requires inner nodes as parameters");

  pll_utree_create_operations(travbuffer,
                              traversal_size,
                              branch_lengths,
                              matrix_indices,
                              operations,
                              &matrix_count,
                              &ops_count);

  pll_update_prob_matrices(partition,
                           params_indices,
                           matrix_indices,
                           branch_lengths,
                           matrix_count);

  pll_update_partials(partition, operations, ops_count);

  logl = pll_compute_edge_loglikelihood(partition,
                                        tree->clv_index,
                                        tree->scaler_index,
                                        tree->back->clv_index,
                                        tree->back->scaler_index,
                                        tree->pmatrix_index,
                                        params_indices,
                                        NULL);

  printf("Recomputed Log-L: %f\n", logl);

  if (fabs(logl - save_logl) < 1e-7)
    printf("Likelihoods OK!\n");
  else
    fatal("Error: Saved and loaded logL do not agree!!\n");

  //pll_utree_show_ascii(tree, (1<<5)-1);

  /* clean */
  pll_partition_destroy(partition);
  pll_utree_graph_destroy(tree, NULL);

  free(block_map);
  free(travbuffer);
  free(branch_lengths);
  free(matrix_indices);
  free(operations);

  return 0;
}
