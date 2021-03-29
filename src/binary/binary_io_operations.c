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

 /**
  * @file bionary_io_operations.c
  *
  * @brief Auxiliary functions for I/O module
  *
  * @author Diego Darriba
  */

#include "binary_io_operations.h"
#include "../pllmod_common.h"

int bin_fwrite(void * data, size_t size, size_t count, FILE * file)
{
  size_t ret = fwrite(data, size, count, file);
  if (ret != count)
  {
    file_io_error(file, ftell(file), "write data");
    return PLL_FAILURE;
  }
  return PLL_SUCCESS;
}

int bin_fread(void * data, size_t size, size_t count, FILE * file)
{
  size_t ret = fread(data, size, count, file);
  if (ret != count)
  {
    file_io_error(file, ftell(file), "read data");
    return PLL_FAILURE;
  }
  return PLL_SUCCESS;
}

int binary_block_header_apply(FILE * bin_file,
                              pll_block_header_t * block_header,
                              int (*bin_func)(void *, size_t, size_t, FILE *))
{
  strcpy(block_header->pad, "000");

  if (!bin_func (block_header, sizeof(pll_block_header_t), 1, bin_file))
  {
    file_io_error(bin_file, ftell(bin_file), "block header apply");
    return PLL_FAILURE;
  }
  return PLL_SUCCESS;
}

int binary_update_header(FILE * bin_file,
                         pll_block_header_t * header)
{
 unsigned int next_block;
 pll_block_map_t next_map;

 long int cur_position = ftell(bin_file);

 pll_binary_header_t bin_header;

 /* update header */
 if (fseek(bin_file, 0, SEEK_SET) == -1)
 {
   file_io_error(bin_file, 0, "update position to head");
   return PLL_FAILURE;
 }

  if (!bin_fread(&bin_header, sizeof(pll_binary_header_t), 1, bin_file))
 {
   file_io_error(bin_file, cur_position, "update binary header [r]");
   return PLL_FAILURE;
 }
 if (fseek(bin_file, 0, SEEK_SET) == -1)
 {
   file_io_error(bin_file, 0, "update position to head");
   return PLL_FAILURE;
 }

 next_block = bin_header.n_blocks;
 ++bin_header.n_blocks;

 if (!bin_fwrite(&bin_header, sizeof(pll_binary_header_t), 1, bin_file))
 {
   file_io_error(bin_file, cur_position, "update binary header [w]");
   return PLL_FAILURE;
 }

  if (header && (header->attributes & PLLMOD_BIN_ATTRIB_UPDATE_MAP))
 {
   /* update map */
   assert(next_block < bin_header.max_blocks);
   if (fseek(bin_file, next_block * sizeof(pll_block_map_t), SEEK_CUR) == -1)
   {
     file_io_error(bin_file, next_block * sizeof(pll_block_map_t),
                   "update position to map");
     return PLL_FAILURE;
   }

   next_map.block_id     = header->block_id;
   next_map.block_offset = cur_position;
   if (!bin_fwrite(&next_map, sizeof(pll_block_map_t), 1, bin_file))
   {
     file_io_error(bin_file, cur_position, "update binary map [w]");
     return PLL_FAILURE;
   }
 }

 /* move back to original position */
 if (fseek(bin_file, cur_position, SEEK_SET) == -1)
 {
   file_io_error(bin_file, cur_position, "update position to block");
   return PLL_FAILURE;
 }

 return PLL_SUCCESS;
}

long int binary_get_offset(FILE *bin_file, int block_id)
{
  pll_block_map_t * map;
  unsigned int i, n_blocks;
  long int offset = PLLMOD_BIN_INVALID_OFFSET;

  map = pllmod_binary_get_map(bin_file, &n_blocks);

  if (!map)
    return PLLMOD_BIN_INVALID_OFFSET;

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

int binary_partition_desc_apply (FILE * bin_file,
                             pll_partition_t * partition,
                             unsigned int attributes,
                             int (*bin_func)(void *, size_t, size_t, FILE *))
{
  PLLMOD_UNUSED(attributes);

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
  bin_func (&partition->asc_bias_alloc, sizeof(int), 1, bin_file);

  /* The variables below are used only if PATTERN_TIP is active. Otherwise
     they could be uninitialized and hence raise a valgrind error if we try
     to write them into the binary file. */
  if (!(partition->attributes & PLL_ATTRIB_PATTERN_TIP))
    partition->maxstates = 0;
  bin_func (&partition->maxstates, sizeof(unsigned int), 1, bin_file);

  return PLL_SUCCESS;
}

int binary_partition_body_apply (FILE * bin_file,
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
  unsigned int sites_alloc = partition->asc_bias_alloc ?
                  partition->sites + partition->states : partition->sites;

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
    if (partition->attributes & PLL_ATTRIB_SITE_REPEATS)
    {
      bin_func (partition->repeats->pernode_ids, sizeof(unsigned int),
          partition->clv_buffers + tips, bin_file);
      bin_func (partition->repeats->perscale_ids, sizeof(unsigned int),
          partition->scale_buffers, bin_file);
    }

  if (attributes & PLLMOD_BIN_ATTRIB_PARTITION_DUMP_CLV)
  {
    unsigned int first_clv_index = 0;

    /* dump tipchars if TIP_PATTERN is used*/
    if (partition->attributes & PLL_ATTRIB_PATTERN_TIP)
    {
      for (i = 0; i < tips; ++i)
      {
        bin_func (partition->tipchars[i],
                  sizeof(unsigned char),
                  sites_alloc,
                  bin_file);
      }
      bin_func (partition->charmap, sizeof(char), PLL_ASCII_SIZE, bin_file);
      first_clv_index = tips;

      unsigned int l2_maxstates = (unsigned int) ceil(log2(
                                                       partition->maxstates));

      /* allocate space for the precomputed tip-tip likelihood vector */
      size_t alloc_size = (1 << (2 * l2_maxstates)) *
                          (partition->states_padded * partition->rate_cats);
      bin_func (partition->ttlookup, sizeof(double), alloc_size, bin_file);
      bin_func (partition->tipmap, sizeof(char), PLL_ASCII_SIZE, bin_file);
    }
    if (partition->attributes & PLL_ATTRIB_SITE_REPEATS) 
    {
      /* dump repeats */
      for (i = first_clv_index; i < (partition->clv_buffers + tips); ++i) 
      {
        unsigned int clvs_to_alloc = pll_get_sites_number(partition, i);
        if (clvs_to_alloc != partition->repeats->pernode_allocated_clvs[i])
        {
          partition->repeats->reallocate_repeats(partition, i,
              (i < tips) ? (unsigned int)PLL_SCALE_BUFFER_NONE : i - tips, clvs_to_alloc); 
        }
        if (partition->repeats->pernode_ids[i]) // repeats
        {
          bin_func(partition->repeats->pernode_site_id[i], sizeof(unsigned int),
            sites_alloc, bin_file);
          bin_func(partition->repeats->pernode_id_site[i], sizeof(unsigned int),
            pll_get_sites_number(partition, i), bin_file);
        } 
      }
    }

    /* dump CLVs and scalers*/
    for (i = first_clv_index; i < (partition->clv_buffers + tips); ++i)
    {
      bin_func (partition->clv[i], sizeof(double),
                pll_get_clv_size(partition, i), bin_file);
    }
    for (i = 0; i < partition->scale_buffers; ++i)
      bin_func (partition->scale_buffer[i], sizeof(unsigned int),
                pll_get_sites_number(partition, partition->tips + i), bin_file);
  }

  if (attributes & PLLMOD_BIN_ATTRIB_PARTITION_DUMP_WGT)
  {
    /* dump pattern weights */
    bin_func (partition->pattern_weights, sizeof(unsigned int), sites_alloc,
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

int binary_partition_apply(FILE * bin_file,
                       pll_partition_t * partition,
                       unsigned int attributes,
                       int (*bin_func)(void *, size_t, size_t, FILE *))
{
  if (!binary_partition_desc_apply(bin_file, partition, attributes, bin_func))
    return PLL_FAILURE;
  if (!binary_partition_body_apply(bin_file, partition, attributes, bin_func))
    return PLL_FAILURE;

  return PLL_SUCCESS;
}

int binary_repeats_apply (FILE * bin_file,
                  pll_partition_t * partition,
                  unsigned int attributes,
                  size_t nodes,
                  int (*bin_func)(void *, size_t, size_t, FILE *))
{
  PLLMOD_UNUSED(attributes);
  if (!bin_func(partition->repeats->pernode_ids, sizeof(unsigned int), nodes, bin_file))
  {
    pllmod_set_error(PLLMOD_BIN_ERROR_LOADSTORE,
                     "Error loading/storing repeats");
    return PLL_FAILURE;
  }
  if (!bin_func(partition->repeats->pernode_allocated_clvs, sizeof(unsigned int), nodes, bin_file))
  {
    pllmod_set_error(PLLMOD_BIN_ERROR_LOADSTORE,
                     "Error loading/storing repeats");
    return PLL_FAILURE;
  }
  return PLL_SUCCESS;
}


int binary_clv_apply (FILE * bin_file,
                      pll_partition_t * partition,
                      unsigned int clv_index,
                      unsigned int attributes,
                      size_t clv_size,
                      int (*bin_func)(void *, size_t, size_t, FILE *))
{
  PLLMOD_UNUSED(attributes);
  if (clv_index > (partition->tips + partition->clv_buffers))
  {
    pllmod_set_error(PLLMOD_BIN_ERROR_INVALID_INDEX,
                     "Invalid CLV index");
    return PLL_FAILURE;
  }
  if (!bin_func(partition->clv[clv_index], sizeof(double), clv_size, bin_file))
  {
    pllmod_set_error(PLLMOD_BIN_ERROR_LOADSTORE,
                     "Error loading/storing CLV");
    return PLL_FAILURE;
  }
  if ((partition->attributes & PLL_ATTRIB_SITE_REPEATS) 
      && partition->repeats->pernode_ids[clv_index]) 
  {
    unsigned int uncompressed_sites = partition->sites + 
      (partition->asc_bias_alloc ? partition->states : 0);
    unsigned int compressed_sites = pll_get_sites_number(partition, clv_index);
    if (!bin_func(partition->repeats->pernode_site_id[clv_index], 
          sizeof(unsigned int), uncompressed_sites, bin_file))
    {
      pllmod_set_error(PLLMOD_BIN_ERROR_LOADSTORE,
                       "Error loading/storing CLV (site_id)");
      return PLL_FAILURE;
    }
    if (!bin_func(partition->repeats->pernode_id_site[clv_index], 
          sizeof(unsigned int), compressed_sites, bin_file))
    {
      pllmod_set_error(PLLMOD_BIN_ERROR_LOADSTORE,
                       "Error loading/storing CLV (id_site)");
      return PLL_FAILURE;
    }
  }
  return PLL_SUCCESS;
}

int binary_node_apply (FILE * bin_file,
                       pll_unode_t * node,
                       int write,
                       int (*bin_func)(void *, size_t, size_t, FILE *))
{
  char * label = 0;
  unsigned long label_len = 0;

  bin_func(node, sizeof(pll_unode_t), 1, bin_file);
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

/* static functions */

void file_io_error (FILE * bin_file, long int setp, const char * msg)
{
  assert(setp >= PLLMOD_BIN_INVALID_OFFSET);

  /* if offset is valid, we apply it */
  if (setp != PLLMOD_BIN_INVALID_OFFSET)
    fseek(bin_file, setp, SEEK_SET);

  /* update error data */
  pllmod_set_error(PLLMOD_BIN_ERROR_LOADSTORE,
                   "Binary file I/O error: %s", msg);
}
