/*
 * binary_io_operations.c
 *
 *  Created on: Mar 23, 2016
 *      Author: diego
 */

#include "binary_io_operations.h"

int bin_fwrite(void * data, size_t size, size_t count, FILE * file)
{
  size_t ret = fwrite(data, size, count, file);
  if (ret != count)
  {
    return PLL_FAILURE;
  }
  return PLL_SUCCESS;
}

int bin_fread(void * data, size_t size, size_t count, FILE * file)
{
  size_t ret = fread(data, size, count, file);
  if (ret != count)
  {
    return PLL_FAILURE;
  }
  return PLL_SUCCESS;
}

long int binary_get_offset(FILE *bin_file, int block_id)
{
  pll_block_map_t * map;
  unsigned int i, n_blocks;
  long int offset = PLL_BINARY_INVALID_OFFSET;

  map = pll_binary_get_map(bin_file, &n_blocks);

  if (!map)
    return PLL_BINARY_INVALID_OFFSET;

  /* search id */
  for (i=0; i<n_blocks; ++i)
  {
    if (map[i].block_id == block_id)
    {
      offset = map[i].block_offset;
      break;
    }
  }

  free(map);
  return offset;
}

int binary_apply_to_partition_desc (FILE * bin_file,
                             pll_partition_t * partition,
                             unsigned int attributes,
                             int (*bin_func)(void *, size_t, size_t, FILE *))
{
  /* partition descriptor */
  bin_func (&partition->tips, sizeof(unsigned int), 1, bin_file);
  bin_func (&partition->clv_buffers, sizeof(unsigned int), 1, bin_file);
  bin_func (&partition->states, sizeof(unsigned int), 1, bin_file);
  bin_func (&partition->sites, sizeof(unsigned int), 1, bin_file);
  bin_func (&partition->rate_matrices, sizeof(unsigned int), 1, bin_file);
  bin_func (&partition->prob_matrices, sizeof(unsigned int), 1, bin_file);
  bin_func (&partition->rate_cats, sizeof(unsigned int), 1, bin_file);
  bin_func (&partition->scale_buffers, sizeof(unsigned int), 1, bin_file);
  bin_func (&partition->attributes, sizeof(unsigned int), 1, bin_file);
  bin_func (&partition->states_padded, sizeof(unsigned int), 1, bin_file);

  /* The variables below are used only if PATTERN_TIP is active. Otherwise
     they could be uninitialized and hence raise a valgrind error if we try
     to write them into the binary file. */
  if (!(partition->attributes & PLL_ATTRIB_PATTERN_TIP))
  {
    partition->maxstates      = 0;
    partition->log2_maxstates = 0;
    partition->log2_states    = 0;
    partition->log2_rates     = 0;
  }

  bin_func (&partition->maxstates, sizeof(unsigned int), 1, bin_file);
  bin_func (&partition->log2_maxstates, sizeof(unsigned int), 1, bin_file);
  bin_func (&partition->log2_states, sizeof(unsigned int), 1, bin_file);
  bin_func (&partition->log2_rates, sizeof(unsigned int), 1, bin_file);

  return PLL_SUCCESS;
}

int binary_apply_to_partition_body (FILE * bin_file,
                             pll_partition_t * partition,
                             unsigned int attributes,
                             int (*bin_func)(void *, size_t, size_t, FILE *))
{
  unsigned int i;
  unsigned int tips = partition->tips;
  unsigned int sites = partition->sites;
  unsigned int rate_cats = partition->rate_cats;
  unsigned int states = partition->states;
  unsigned int states_padded = partition->states_padded;
  unsigned int n_subst_rates = states * (states - 1) / 2;
  unsigned int prob_matrices = partition->prob_matrices;
  unsigned int rate_matrices = partition->rate_matrices;

  bin_func (partition->eigen_decomp_valid, sizeof(int), rate_matrices,
            bin_file);
  for (i = 0; i < rate_matrices; ++i)
    bin_func (partition->eigenvecs[i], sizeof(double),
              states * states_padded, bin_file);
  for (i = 0; i < rate_matrices; ++i)
    bin_func (partition->inv_eigenvecs[i], sizeof(double),
              states * states_padded, bin_file);
  for (i = 0; i < rate_matrices; ++i)
    bin_func (partition->eigenvals[i], sizeof(double), states_padded,
              bin_file);
  bin_func (partition->pmatrix[0],
            sizeof(double),
            prob_matrices * states * states_padded * rate_cats, bin_file);
  for (i = 0; i < rate_matrices; ++i)
    bin_func (partition->subst_params[i], sizeof(double),
              n_subst_rates, bin_file);
  for (i = 0; i < rate_matrices; ++i)
    bin_func (partition->frequencies[i], sizeof(double), states_padded,
              bin_file);
  bin_func (partition->rates, sizeof(double), rate_cats, bin_file);
  bin_func (partition->rate_weights, sizeof(double), rate_cats, bin_file);
  bin_func (partition->prop_invar, sizeof(double), rate_matrices, bin_file);

  if (attributes & PLL_BINARY_ATTRIB_PARTITION_DUMP_CLV)
  {
    unsigned int first_clv_index = 0;
    /* dump tipchars if TIP_PATTERN is used*/
    if (partition->attributes & PLL_ATTRIB_PATTERN_TIP)
    {
      for (i = 0; i < tips; ++i)
        bin_func (partition->tipchars[i], sizeof(char), sites, bin_file);
      bin_func (partition->charmap, sizeof(char), PLL_ASCII_SIZE, bin_file);
      first_clv_index = tips;
    }
    /* dump CLVs and scalers*/
    for (i = first_clv_index; i < (partition->clv_buffers + tips); ++i)
    {
      bin_func (partition->clv[i], sizeof(double),
                states_padded * rate_cats * sites, bin_file);
    }
    for (i = 0; i < partition->scale_buffers; ++i)
      bin_func (partition->scale_buffer[i], sizeof(unsigned int), sites,
                bin_file);
  }

  if (attributes & PLL_BINARY_ATTRIB_PARTITION_DUMP_WGT)
  {
    /* dump pattern weights */
    bin_func (partition->pattern_weights, sizeof(unsigned int), sites,
              bin_file);
  }

  for (i = 0; i < rate_matrices; ++i)
    if (partition->prop_invar[i] > 0)
    {
      if (!partition->invariant)
        partition->invariant = (int *)malloc(partition->sites * sizeof(int));
      bin_func (partition->invariant, sizeof(int), sites, bin_file);
      break;
    }

  return PLL_SUCCESS;
}

int binary_apply_to_partition(FILE * bin_file,
                       pll_partition_t * partition,
                       unsigned int attributes,
                       int (*bin_func)(void *, size_t, size_t, FILE *))
{
  if (!binary_apply_to_partition_desc(bin_file, partition, attributes, bin_func))
    return PLL_FAILURE;

  if (!binary_apply_to_partition_body(bin_file, partition, attributes, bin_func))
    return PLL_FAILURE;

  return PLL_SUCCESS;
}

int binary_apply_to_clv (FILE * bin_file,
                  pll_partition_t * partition,
                  unsigned int clv_index,
                  unsigned int attributes,
                  size_t clv_size,
                  int (*bin_func)(void *, size_t, size_t, FILE *))
{
  if (clv_index < 0 || clv_index > (partition->tips + partition->clv_buffers))
  {
    pll_errno = PLL_ERROR_INVALID_INDEX;
    snprintf(pll_errmsg, 200, "Invalid CLV index");
    return PLL_FAILURE;
  }

  if (!bin_func(partition->clv[clv_index], sizeof(double), clv_size, bin_file))
  {
    pll_errno = PLL_ERROR_LOADSTORE;
    snprintf(pll_errmsg, 200, "Error loading/storing CLV");
    return PLL_FAILURE;
  }

  return PLL_SUCCESS;
}

int binary_apply_to_node (FILE * bin_file,
                          pll_utree_t * node,
                          int write,
                          int (*bin_func)(void *, size_t, size_t, FILE *))
{
  char * label = 0;
  unsigned long label_len = 0;

  bin_func(node, sizeof(pll_utree_t), 1, bin_file);
  if (write && node->label)
    label_len = strlen(node->label);
  bin_func(&label_len, sizeof(unsigned long), 1, bin_file);
  if (label_len)
  {
    if (write)
    {
      label = node->label;
    }
    else
    {
      label = (char *) malloc(label_len+1);
      node->label = label;
    }
    bin_func(label, sizeof(char), label_len, bin_file);
    node->label[label_len] = '\0';
  }

  return PLL_SUCCESS;
}
