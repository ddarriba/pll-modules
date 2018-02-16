/*
 Copyright (C) 2017 Diego Darriba, Pierre Barbera

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
 * @file pll_binary.c
 *
 * @brief Binary I/O operations for PLL
 *
 * @author Diego Darriba
 * @author Pierre Barbera
 */

#include "pll_binary.h"
#include "binary_io_operations.h"
#include "../pllmod_common.h"
#include <stdio.h>

static unsigned int get_current_alignment( unsigned int attributes );
static int cb_full_traversal(pll_unode_t * node);

/**
 *  Open file for writing
 *
 *  @param[in] filename file to write to
 *  @param[out] header file header
 *  @param access_type PLLMOD_BIN_ACCESS_[SEQUENTIAL|RANDOM]
 *  @param n_blocks actual or maximum number of blocks if access is random
 *
 *  @return pointer to the file
 */
PLL_EXPORT FILE * pllmod_binary_create(const char * filename,
                                       pll_binary_header_t * header,
                                       unsigned int access_type,
                                       unsigned int n_blocks)
{
  FILE * file = NULL;

  memset(header, 0, sizeof(pll_binary_header_t));
  header->access_type = access_type;
  header->max_blocks = n_blocks;
  header->map_offset = (access_type == PLLMOD_BIN_ACCESS_RANDOM)?
      n_blocks * sizeof(pll_block_map_t):0;
  header->n_blocks = 0;

  if (access_type == PLLMOD_BIN_ACCESS_RANDOM && n_blocks <= 0)
  {
    pllmod_set_error(PLLMOD_BIN_ERROR_INVALID_SIZE,
             "Number of blocks for random access must be greater than 0");
    return NULL;
  }

  file = fopen(filename, "w+b");

  if (!file)
  {
    pllmod_set_error(PLL_ERROR_FILE_OPEN, "Cannot open file for writing");
    return NULL;
  }

  if (!bin_fwrite(header, sizeof(pll_binary_header_t), 1, file))
  {
    pllmod_set_error(PLLMOD_BIN_ERROR_BINARY_IO,
                     "Error writing header to file");
    fclose(file);
    return NULL;
  }

  if(fseek(file, header->map_offset, SEEK_CUR) == -1)
  {
    pllmod_set_error(PLLMOD_BIN_ERROR_BINARY_IO,
                     "Error seeking through file during creation");
    fclose(file);
    return NULL;
  }

  return file;
}

/**
 *  Open file for reading
 *
 *  @param[in] filename file to read from
 *  @param[out] header file header
 *
 *  @return pointer to the file
 */
PLL_EXPORT FILE * pllmod_binary_open(const char * filename,
                                     pll_binary_header_t * header)
{
  FILE * file;

  file = fopen(filename, "rb");

  if (!file)
  {
    pllmod_set_error(PLL_ERROR_FILE_OPEN, "Cannot open file for reading");
    return NULL;
  }

  if (!bin_fread(header, sizeof(pll_binary_header_t), 1, file))
  {
    pllmod_set_error(PLLMOD_BIN_ERROR_BINARY_IO,
                     "Error reading header from file");
    fclose(file);
    return NULL;
  }

  fseek(file, header->map_offset, SEEK_CUR);

  return file;
}

/**
 *  Open file for appending
 *
 *  @param[in] filename file to append to
 *  @param[out] header file header
 *
 *  @return pointer to the file
 */
FILE * pllmod_binary_append_open(const char * filename,
                                 pll_binary_header_t * header)
{
  FILE * file;
  long fpos;

  file = fopen(filename, "r+b");

  if (!file)
  {
    pllmod_set_error(PLL_ERROR_FILE_OPEN, "Cannot open file for appending");
    return NULL;
  }

  if (!bin_fread(header, sizeof(pll_binary_header_t), 1, file))
  {
    pllmod_set_error(PLLMOD_BIN_ERROR_BINARY_IO,
                     "Error reading header from file");
    fclose(file);
    return NULL;
  }

  if (fseek(file, 0, SEEK_END) == -1)
  {
    file_io_error(file, ftell(file), "update position to EOF");
    fclose(file);
    return NULL;
  }
  fpos = ftell(file);
  if (header->map_offset > fpos)
    fpos = header->map_offset;
  if (fseek(file, fpos, SEEK_SET) == -1)
  {
    file_io_error(file, fpos, "update position to last block");
    fclose(file);
    return NULL;
  }

  return file;
}

/**
 *  Closes the binary file
 *
 *  @param[in] bin_file the file
 *
 *  @return true, if OK
 */
PLL_EXPORT int pllmod_binary_close(FILE * bin_file)
{
  return fclose(bin_file);
}

/**
 *  Save a partition to the binary file
 *
 *  @param[in] bin_file binary file
 *  @param[in] block_id id of the block for random access, or local id
 *  @param[in] partition the saved partition
 *  @param[in] attributes the dumped attributes
 *
 * @return PLL_SUCCESS if the data was correctly saved
 *         PLL_FAILURE otherwise (check pll_errmsg for details)
 */

PLL_EXPORT int pllmod_binary_partition_dump(FILE * bin_file,
                                            int block_id,
                                            pll_partition_t * partition,
                                            unsigned int attributes)
{
  pll_block_header_t block_header;
  // unsigned long partition_len = partition_size(partition),
  //               clv_len = 0,
  //               wgt_len = 0;
  long int start_pos = ftell(bin_file),
           end_pos;

  /* fill block header */
  block_header.block_id   = block_id;
  block_header.type       = PLLMOD_BIN_BLOCK_PARTITION;
  block_header.attributes = attributes;
  block_header.block_len  = 0; //partition_len;
  block_header.alignment  = 0;

  /* update main header */
  if (!binary_update_header(bin_file, &block_header))
  {
    return PLL_FAILURE;
  }

  /* dump header */
  if (!binary_block_header_apply(bin_file, &block_header, &bin_fwrite))
  {
    return PLL_FAILURE;
  }

  /* dump data */
  if (!binary_partition_apply(bin_file, partition, attributes, &bin_fwrite))
  {
    return PLL_FAILURE;
  }

  end_pos = ftell(bin_file);

  /* update header */
  if (fseek(bin_file, start_pos, SEEK_SET) == -1)
  {
    file_io_error(bin_file, start_pos, "update position to header");
    return PLL_FAILURE;
  }

  block_header.block_len = (size_t) (end_pos - start_pos);
  if (!binary_block_header_apply(bin_file, &block_header, &bin_fwrite))
  {
    assert(pll_errno);
    return PLL_FAILURE;
  }

  if (fseek(bin_file, end_pos, SEEK_SET) == -1)
  {
    file_io_error(bin_file, end_pos, "update position to end");
    return PLL_FAILURE;
  }

  return PLL_SUCCESS;
}

/**
 *  Load a partition from the binary file
 *
 *  @param[in] bin_file binary file
 *  @param[in] block_id id of the block for random access
 *  @param[in,out] partition if NULL, creates a new partition
 *  @param[in,out] attributes the loaded attributes. If
 *                 PLLMOD_BIN_ATTRIB_PARTITION_LOAD_SKELETON passed here,
 *                 only pointers to the CLVs, tipchars and scalers are
 *                 allocated instead ofthe full memory. Pointers will
 *                 be initialized with NULL
 *  @param offset offset to the data block, if known
 *                0, if access is sequential
 *                PLLMOD_BIN_ACCESS_SEEK, for searching in the file header
 *
 *  @return pointer to the updated (or new) partition
 */

PLL_EXPORT pll_partition_t * pllmod_binary_partition_load(FILE * bin_file,
                                                          int block_id,
                                                          pll_partition_t * partition,
                                                          unsigned int * attributes,
                                                          long int offset)
{
  pll_block_header_t block_header;
  pll_partition_t * local_partition;
  assert(offset >= 0 || offset == PLLMOD_BIN_ACCESS_SEEK);
  unsigned int sites_alloc;
  unsigned int i;

  const int load_skeleton =
    *attributes & PLLMOD_BIN_ATTRIB_PARTITION_LOAD_SKELETON;

  if (offset != 0)
  {
    if (offset == PLLMOD_BIN_ACCESS_SEEK)
    {
      /* find offset */
      offset = binary_get_offset (bin_file, block_id);
      if (offset == PLLMOD_BIN_INVALID_OFFSET)
      {
        pllmod_set_error(PLLMOD_BIN_ERROR_BINARY_IO,
                       "Cannot retrieve offset for block %d", block_id);
        return NULL;
      }
    }

    /* apply offset */
    fseek (bin_file, offset, SEEK_SET);
  }

  if (!binary_block_header_apply(bin_file, &block_header, &bin_fread))
    return NULL;

  if (block_header.type != PLLMOD_BIN_BLOCK_PARTITION)
  {
    pllmod_set_error(PLLMOD_BIN_ERROR_BLOCK_MISMATCH,
                     "Block type is %d and should be %d",
                     block_header.type, PLLMOD_BIN_BLOCK_PARTITION);
    return NULL;
  }

  *attributes = block_header.attributes;

  if (partition)
  {
    local_partition = partition;
  }
  else
  {
    /* create new */
    pll_partition_t aux_partition;
    if (!binary_partition_desc_apply (bin_file,
                                      &aux_partition,
                                      *attributes,
                                      &bin_fread))
    {
      return NULL;
    }

    unsigned int clv_buffers = load_skeleton ? 1 : aux_partition.clv_buffers;
    unsigned int tips = load_skeleton ? 0 : aux_partition.tips;
    unsigned int scale_buffers = load_skeleton ? 1 : aux_partition.scale_buffers;
    local_partition = pll_partition_create(
        tips,
        clv_buffers,
        aux_partition.states,
        aux_partition.sites,
        aux_partition.rate_matrices,
        aux_partition.prob_matrices,
        aux_partition.rate_cats,
        scale_buffers,
        aux_partition.attributes);

    if (load_skeleton)
    {
      if (local_partition->clv)
      {
        size_t start = (local_partition->attributes & PLL_ATTRIB_PATTERN_TIP) ?
                        local_partition->tips : 0;
        for (i = start; i < local_partition->clv_buffers + local_partition->tips; ++i)
          pll_aligned_free(local_partition->clv[i]);
      }
      free(local_partition->clv);
      local_partition->clv_buffers = aux_partition.clv_buffers;
      local_partition->tips = aux_partition.tips;
      local_partition->nodes = aux_partition.clv_buffers + aux_partition.tips;;
      local_partition->clv = (double**) calloc(local_partition->clv_buffers +
                                               local_partition->tips,
                                               sizeof(double*));
      if (local_partition->scale_buffer)
        for (i = 0; i < local_partition->scale_buffers; ++i)
      free(local_partition->scale_buffer[i]);
      free(local_partition->scale_buffer);
      local_partition->scale_buffers = aux_partition.scale_buffers;
      local_partition->scale_buffer = (unsigned int **) calloc(
                                                local_partition->scale_buffers,
                                                sizeof(unsigned int *));

      // manually set the tips so that the rest of the code callocs correctly
      local_partition->tips = aux_partition.tips;
      
      if (local_partition->attributes & PLL_ATTRIB_SITE_REPEATS) 
      {
        free(local_partition->repeats->pernode_ids);
        free(local_partition->repeats->pernode_allocated_clvs);
        free(local_partition->repeats->perscale_ids);
        free(local_partition->repeats->pernode_site_id[0]);
        free(local_partition->repeats->pernode_id_site[0]);
        free(local_partition->repeats->pernode_site_id);
        free(local_partition->repeats->pernode_id_site);
        local_partition->repeats->pernode_ids = calloc(local_partition->clv_buffers 
            + local_partition->tips, sizeof(unsigned int));
        local_partition->repeats->pernode_allocated_clvs = calloc(local_partition->clv_buffers 
            + local_partition->tips, sizeof(unsigned int));
        local_partition->repeats->perscale_ids = calloc(local_partition->scale_buffers,
            sizeof(unsigned int));
        local_partition->repeats->pernode_site_id = calloc(local_partition->clv_buffers 
            + local_partition->tips, sizeof(unsigned int *));
        local_partition->repeats->pernode_id_site = calloc(local_partition->clv_buffers 
            + local_partition->tips, sizeof(unsigned int *));
      }
    }

    /* initialize extra variables */
    local_partition->maxstates = aux_partition.maxstates;
    local_partition->asc_bias_alloc = aux_partition.asc_bias_alloc;

    sites_alloc = local_partition->asc_bias_alloc ?
                   local_partition->sites + local_partition->states :
                   local_partition->sites;
    if (local_partition->attributes & PLL_ATTRIB_PATTERN_TIP)
    {
      /* allocate tip character arrays */
      local_partition->tipchars =
        (unsigned char **)calloc(local_partition->tips,
                                 sizeof(unsigned char *));
      if (!local_partition->tipchars)
      {
        pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                  "Cannot allocate space for storing tip characters.");
        return PLL_FAILURE;
      }

      if (!(local_partition->charmap = (unsigned char *)calloc(PLL_ASCII_SIZE,
                                                        sizeof(unsigned char))))
      {
        pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                  "Cannot allocate charmap for tip-tip precomputation.");
        return PLL_FAILURE;
      }

      if (!(local_partition->tipmap = (pll_state_t *)calloc(PLL_ASCII_SIZE,
                                                       sizeof(pll_state_t))))
      {
        pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                  "Cannot allocate tipmap for tip-tip precomputation.");
        return PLL_FAILURE;
      }

      if (!load_skeleton)
      {
        for (i = 0; i < local_partition->tips ; ++i)
        {
          local_partition->tipchars[i] = (unsigned char *)malloc(sites_alloc *
                                                           sizeof(unsigned char));
          if (!local_partition->tipchars[i])
          {
            pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                      "Cannot allocate space for storing tip characters.");
            return PLL_FAILURE;
          }
        }
      }

      if ((local_partition->states == 4) &&
         (local_partition->attributes & PLL_ATTRIB_ARCH_AVX))
      {
        local_partition->ttlookup = pll_aligned_alloc(1024 *
                                                local_partition->rate_cats *
                                                sizeof(double),
                                                local_partition->alignment);
      }
      else
      {
        unsigned int l2_maxstates =
          (unsigned int) ceil(log2(local_partition->maxstates));
        size_t alloc_size = (1 << (2 * l2_maxstates)) *
                            (local_partition->states_padded *
                            local_partition->rate_cats);
        local_partition->ttlookup = pll_aligned_alloc(alloc_size *
                                                      sizeof(double),
                                                    local_partition->alignment);
      }
      if (!local_partition->ttlookup)
      {
        pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                "Cannot allocate space for storing precomputed tip-tip CLVs.");
        return PLL_FAILURE;
      }
    }

    if (!local_partition)
    {
      return NULL;
    }
  }

  if (!binary_partition_body_apply (bin_file,
                                    local_partition,
                                    *attributes,
                                    &bin_fread))
  {
    pll_partition_destroy(local_partition);
    return NULL;
  }

  return local_partition;
}

/**
 *  Save repeats to the binary file
 *
 *  @param[in] bin_file binary file
 *  @param[in] block_id id of the block for random access, or local id
 *  @param[in] partition the partition containing the saved CLV
 *  @param[in] attributes the dumped attributes
 *
 *  @return PLL_SUCCESS if the data was correctly saved
 *          PLL_FAILURE otherwise (check pll_errmsg for details)
 */
PLL_EXPORT int pllmod_binary_repeats_dump(FILE * bin_file,
                                      int block_id,
                                      pll_partition_t * partition,
                                      unsigned int attributes)
{
  assert(partition);
  if (!(partition->attributes & PLL_ATTRIB_SITE_REPEATS)) 
  {
    return PLL_SUCCESS;
  }
  assert(partition->repeats);
  int retval;
  pll_block_header_t block_header;

  size_t nodes = partition->tips + partition->clv_buffers;
  /* fill block header */
  block_header.block_id   = block_id;
  block_header.type       = PLLMOD_BIN_BLOCK_REPEATS;
  block_header.attributes = attributes;
  block_header.block_len  = nodes * sizeof(unsigned int);
  block_header.alignment  = 0;

  /* update main header */
  if(!binary_update_header(bin_file, &block_header))
  {
    return PLL_FAILURE;
  }

  /* dump block header */
  if (!binary_block_header_apply(bin_file, &block_header, &bin_fwrite))
    return PLL_FAILURE;

  /* dump data */
  retval = binary_repeats_apply (bin_file,
                             partition,
                             attributes,
                             nodes,
                             &bin_fwrite);
  return retval;
}

/**
 *  Load a CLV from the binary file
 *
 *  @param[in] bin_file binary file
 *  @param[in] block_id id of the block for random access
 *  @param[in,out] partition the partition where the CLV will be stored
 *  @param[out] attributes the loaded attributes
 *  @param offset offset to the data block, if known
 *                0, if access is sequential
 *                PLLMOD_BIN_ACCESS_SEEK, for searching in the file header
 *
 *  @return PLL_SUCCESS if the data was correctly loaded
 *          PLL_FAILURE otherwise (check pll_errmsg for details)
 */

PLL_EXPORT int pllmod_binary_repeats_load(FILE * bin_file,
                                      int block_id,
                                      pll_partition_t * partition,
                                      unsigned int * attributes,
                                      long int offset)
{
  if (!(partition->attributes & PLL_ATTRIB_SITE_REPEATS)) 
  {
    return PLL_SUCCESS;
  }
  int retval;
  pll_block_header_t block_header;
  assert(partition);
  assert(offset >= 0 || offset == PLLMOD_BIN_ACCESS_SEEK);
  assert(partition->repeats);

  if (offset != 0)
  {
    if (offset == PLLMOD_BIN_ACCESS_SEEK)
    {
      /* find offset */
      offset = binary_get_offset (bin_file, block_id);
      if (offset == PLLMOD_BIN_INVALID_OFFSET)
      {
        pllmod_set_error(PLLMOD_BIN_ERROR_MISSING_BLOCK,
                      "Cannot retrieve offset for block %d", block_id);
        return PLL_FAILURE;
      }
    }
    /* apply offset */
    fseek (bin_file, offset, SEEK_SET);
  }

  /* read and validate header */
  if (!binary_block_header_apply(bin_file, &block_header, &bin_fread))
    return PLL_FAILURE;

  if (block_header.type != PLLMOD_BIN_BLOCK_REPEATS)
  {
    pllmod_set_error(PLLMOD_BIN_ERROR_BLOCK_MISMATCH,
                  "Block type is %d and should be %d",
                  block_header.type, PLLMOD_BIN_BLOCK_REPEATS);
    return PLL_FAILURE;
  }

  size_t nodes = partition->tips + partition->clv_buffers;
  if (block_header.block_len != (nodes * sizeof(unsigned int)))
  {
      pllmod_set_error(PLLMOD_BIN_ERROR_BLOCK_LENGTH,
                    "Wrong block length");
      return PLL_FAILURE;
  }
  *attributes = block_header.attributes;
  retval = binary_repeats_apply (bin_file,
                         partition,
                         *attributes,
                         nodes,
                         &bin_fread);
  return retval;
}


/**
 *  Save a repeat entry to the binary file
 *
 *  @param[in] bin_file binary file
 *  @param[in] block_id id of the block for random access, or local id
 *  @param[in] partition the partition containing the saved repeats
 *  @param[in] clv_index the index of the repeat entry
 *  @param[in] attributes the dumped attributes
 *
 *  @return PLL_SUCCESS if the data was correctly saved
 *          PLL_FAILURE otherwise (check pll_errmsg for details)
 */
PLL_EXPORT int pllmod_binary_pernoderepeats_dump(FILE * bin_file,
                                      int block_id,
                                      pll_partition_t * partition,
                                      unsigned int clv_index,
                                      unsigned int attributes)
{
  int retval;
  pll_block_header_t block_header;
  unsigned int sites = partition->attributes & PLL_ATTRIB_SITE_REPEATS ?
    partition->repeats->pernode_ids[clv_index] : partition->sites;
  unsigned int sites_alloc = partition->asc_bias_alloc ?
                 sites + partition->states :
                 sites;

  size_t clv_size = sites_alloc * partition->states_padded *
                      partition->rate_cats;
  /* fill block header */
  block_header.block_id   = block_id;
  block_header.type       = PLLMOD_BIN_BLOCK_CLV;
  block_header.attributes = attributes;
  block_header.block_len  = clv_size * sizeof(double);
  block_header.alignment  = 0;

  /* update main header */
  if(!binary_update_header(bin_file, &block_header))
  {
    return PLL_FAILURE;
  }

  /* dump block header */
  if (!binary_block_header_apply(bin_file, &block_header, &bin_fwrite))
    return PLL_FAILURE;

  /* dump data */
  retval = binary_clv_apply (bin_file,
                             partition,
                             clv_index,
                             attributes,
                             clv_size,
                             &bin_fwrite);

  return retval;
}

/**
 *  Load a repeat entry from the binary file
 *
 *  @param[in] bin_file binary file
 *  @param[in] block_id id of the block for random access
 *  @param[in,out] partition the partition where the repeats will be stored
 *  @param[in] clv_index index of the repeat entry
 *  @param[out] attributes the loaded attributes
 *  @param offset offset to the data block, if known
 *                0, if access is sequential
 *                PLLMOD_BIN_ACCESS_SEEK, for searching in the file header
 *
 *  @return PLL_SUCCESS if the data was correctly loaded
 *          PLL_FAILURE otherwise (check pll_errmsg for details)
 */

PLL_EXPORT int pllmod_binary_pernoderepeats_load(FILE * bin_file,
                                      int block_id,
                                      pll_partition_t * partition,
                                      unsigned int clv_index,
                                      unsigned int * attributes,
                                      long int offset)
{
  
  return 0;
}


/**
 *  Save a CLV to the binary file
 *
 *  @param[in] bin_file binary file
 *  @param[in] block_id id of the block for random access, or local id
 *  @param[in] partition the partition containing the saved CLV
 *  @param[in] clv_index the index of the CLV
 *  @param[in] attributes the dumped attributes
 *
 *  @return PLL_SUCCESS if the data was correctly saved
 *          PLL_FAILURE otherwise (check pll_errmsg for details)
 */
PLL_EXPORT int pllmod_binary_clv_dump(FILE * bin_file,
                                      int block_id,
                                      pll_partition_t * partition,
                                      unsigned int clv_index,
                                      unsigned int attributes)
{
  int retval;
  pll_block_header_t block_header;
  size_t clv_size = pll_get_clv_size(partition, clv_index);
  /* fill block header */
  block_header.block_id   = block_id;
  block_header.type       = PLLMOD_BIN_BLOCK_CLV;
  block_header.attributes = attributes;
  block_header.block_len  = clv_size * sizeof(double);
  block_header.alignment  = 0;
  
  if ((partition->attributes & PLL_ATTRIB_SITE_REPEATS) 
      && partition->repeats->pernode_ids[clv_index]) 
  {
    unsigned int uncompressed_sites = partition->sites + 
      (partition->asc_bias_alloc ? partition->states : 0);
    unsigned int compressed_sites = pll_get_sites_number(partition, clv_index);
    block_header.block_len += (uncompressed_sites + compressed_sites) * sizeof(unsigned int);
  }
  /* update main header */
  if(!binary_update_header(bin_file, &block_header))
  {
    return PLL_FAILURE;
  }

  /* dump block header */
  if (!binary_block_header_apply(bin_file, &block_header, &bin_fwrite))
    return PLL_FAILURE;

  /* dump data */
  retval = binary_clv_apply (bin_file,
                             partition,
                             clv_index,
                             attributes,
                             clv_size,
                             &bin_fwrite);

  return retval;
}




/**
 *  Load a CLV from the binary file
 *
 *  @param[in] bin_file binary file
 *  @param[in] block_id id of the block for random access
 *  @param[in,out] partition the partition where the CLV will be stored
 *  @param[in] clv_index index of the CLV
 *  @param[out] attributes the loaded attributes
 *  @param offset offset to the data block, if known
 *                0, if access is sequential
 *                PLLMOD_BIN_ACCESS_SEEK, for searching in the file header
 *
 *  @return PLL_SUCCESS if the data was correctly loaded
 *          PLL_FAILURE otherwise (check pll_errmsg for details)
 */

PLL_EXPORT int pllmod_binary_clv_load(FILE * bin_file,
                                      int block_id,
                                      pll_partition_t * partition,
                                      unsigned int clv_index,
                                      unsigned int * attributes,
                                      long int offset)
{
  int retval;
  pll_block_header_t block_header;

  assert (partition);
  assert(offset >= 0 || offset == PLLMOD_BIN_ACCESS_SEEK);
  
  size_t clv_size = pll_get_clv_size(partition, clv_index);
  if (offset != 0)
  {
    if (offset == PLLMOD_BIN_ACCESS_SEEK)
    {
      /* find offset */
      offset = binary_get_offset (bin_file, block_id);
      if (offset == PLLMOD_BIN_INVALID_OFFSET)
      {
        pllmod_set_error(PLLMOD_BIN_ERROR_MISSING_BLOCK,
                      "Cannot retrieve offset for block %d", block_id);
        return PLL_FAILURE;
      }
    }

    /* apply offset */
    fseek (bin_file, offset, SEEK_SET);
  }

  /* read and validate header */
  if (!binary_block_header_apply(bin_file, &block_header, &bin_fread))
    return PLL_FAILURE;

  if (block_header.type != PLLMOD_BIN_BLOCK_CLV)
  {
    pllmod_set_error(PLLMOD_BIN_ERROR_BLOCK_MISMATCH,
                  "Block type is %d and should be %d",
                  block_header.type, PLLMOD_BIN_BLOCK_CLV);
    return PLL_FAILURE;
  }
  
  unsigned int block_len = clv_size * sizeof(double);

  if ((partition->attributes & PLL_ATTRIB_SITE_REPEATS) 
      && partition->repeats->pernode_ids[clv_index]) 
  {
    unsigned int uncompressed_sites = partition->sites + 
      (partition->asc_bias_alloc ? partition->states : 0);
    unsigned int compressed_sites = pll_get_sites_number(partition, clv_index);
    block_len += (uncompressed_sites + compressed_sites) * sizeof(unsigned int);
    free(partition->repeats->pernode_site_id[clv_index]);
    free(partition->repeats->pernode_id_site[clv_index]);
    partition->repeats->pernode_site_id[clv_index] = malloc(uncompressed_sites * sizeof(unsigned int));
    partition->repeats->pernode_id_site[clv_index] = malloc(compressed_sites * sizeof(unsigned int));
  }

  if (block_header.block_len != block_len)
  {
      pllmod_set_error(PLLMOD_BIN_ERROR_BLOCK_LENGTH,
                    "Wrong block length");
      return PLL_FAILURE;
  }

  *attributes = block_header.attributes;

  retval = binary_clv_apply (bin_file,
                         partition,
                         clv_index,
                         *attributes,
                         clv_size,
                         &bin_fread);

  return retval;
}

/**
 *  Save an unrooted tree to the binary file
 *
 *  @param[in] bin_file binary file
 *  @param[in] block_id id of the block for random access, or local id
 *  @param[in] tree the unrooted tree to be dumped
 *  @param[in] tip_count the number of tips in the tree
 *  @param[in] attributes the loaded attributes
 *
 *  @return PLL_SUCCESS if the data was correctly saved
 *          PLL_FAILURE otherwise (check pll_errmsg for details)
 */
PLL_EXPORT int pllmod_binary_utree_dump(FILE * bin_file,
                                        int block_id,
                                        pll_unode_t * tree,
                                        unsigned int tip_count,
                                        unsigned int attributes)
{
  pll_unode_t ** travbuffer;
  pll_block_header_t block_header;
  unsigned int i, n_nodes, n_inner, n_utrees, trav_size;
  int retval;

  /* reset error */
  pll_errno = 0;

  n_inner = tip_count - 2;
  n_nodes = tip_count + n_inner;
  n_utrees = tip_count + 3 * n_inner;

  travbuffer = (pll_unode_t **)malloc(n_nodes* sizeof(pll_unode_t *));

  if (!tree->next)
  {
    pllmod_set_error(PLLMOD_BIN_ERROR_BINARY_IO,
                     "Tree should not be a tip node");
    return PLL_FAILURE;
  }
  if (!pll_utree_traverse(tree,
                          PLL_TREE_TRAVERSE_POSTORDER,
                          cb_full_traversal,
                          travbuffer,
                          &trav_size))
  {
    pllmod_set_error(PLLMOD_BIN_ERROR_BINARY_IO,
                     "Error traversing utree");
    return PLL_FAILURE;
  }

  assert (trav_size == n_nodes);

  block_header.block_id   = block_id;
  block_header.type       = PLLMOD_BIN_BLOCK_TREE;
  block_header.attributes = attributes;
  block_header.block_len  = n_utrees * sizeof(pll_unode_t);
  block_header.alignment  = 0;

  /* update main header */
  if(!binary_update_header(bin_file, &block_header))
  {
    assert(pll_errno);
    return PLL_FAILURE;
  }

  /* dump block header */
  if (!binary_block_header_apply(bin_file, &block_header, &bin_fwrite))
  {
    assert(pll_errno);
    return PLL_FAILURE;
  }

  /* traverse and dump data */
  for (i=0; i<trav_size; ++i)
  {
    if (!binary_node_apply (bin_file,
                               travbuffer[i],
                               1,
                               bin_fwrite))
    {
      assert(pll_errno);
      return PLL_FAILURE;
    }

    if (travbuffer[i]->next)
    {
      retval = binary_node_apply (bin_file,
                                travbuffer[i]->next,
                                1,
                                bin_fwrite);
      retval &= binary_node_apply (bin_file,
                                 travbuffer[i]->next->next,
                                 1,
                                 bin_fwrite);
      if (!retval)
      {
        assert(pll_errno);
        return PLL_FAILURE;
      }
    }
  }

  free(travbuffer);

  return PLL_SUCCESS;
}

/**
 *  Load an unrooted tree from the binary file
 *
 *  @param[in] bin_file binary file
 *  @param[in] block_id id of the block for random access
 *  @param[out] attributes the block attributes
 *  @param offset offset to the data block, if known
 *                0, if access is sequential
 *                PLLMOD_BIN_ACCESS_SEEK, for searching in the file header
 *
 *  @return pointer to the updated (or new) partition
 */
PLL_EXPORT pll_unode_t * pllmod_binary_utree_load(FILE * bin_file,
                                                  int block_id,
                                                  unsigned int * attributes,
                                                  long int offset)
{
  unsigned int i, n_tips, n_tip_check, n_nodes;
  long n_utrees;
  pll_block_header_t block_header;
  pll_unode_t ** tree_stack;
  pll_unode_t * tree;
  unsigned int tree_stack_top;
  int retval;

  assert(offset >= 0 || offset == PLLMOD_BIN_ACCESS_SEEK);

  if (offset != 0)
  {
    if (offset == PLLMOD_BIN_ACCESS_SEEK)
    {
      /* find offset */
      offset = binary_get_offset (bin_file, block_id);
      if (offset == PLLMOD_BIN_INVALID_OFFSET)
      {
        pllmod_set_error(PLLMOD_BIN_ERROR_BINARY_IO,
                         "Cannot retrieve offset for block %d", block_id);
        return NULL;
      }
    }

    /* apply offset */
    fseek (bin_file, offset, SEEK_SET);
  }

  /* read and validate header */
  if (!binary_block_header_apply(bin_file, &block_header, &bin_fread))
  {
    assert(pll_errno);
    return NULL;
  }

  if (block_header.type != PLLMOD_BIN_BLOCK_TREE)
  {
    pllmod_set_error(PLLMOD_BIN_ERROR_BLOCK_MISMATCH,
                     "Block type is %d and should be %d",
                     block_header.type, PLLMOD_BIN_BLOCK_TREE);
    return NULL;
  }

  n_utrees = block_header.block_len/sizeof(pll_unode_t);
  n_tips   = (unsigned int) ((n_utrees + 6) / 4);
  n_nodes  = 2*n_tips - 2;
  assert( n_utrees % 4 == 2 );

  *attributes = block_header.attributes;

  /* allocate stack for at most 'n_tips' nodes */
  tree_stack = (pll_unode_t **) malloc(n_tips * sizeof (pll_unode_t *));
  tree_stack_top = 0;

  /* read nodes */
  n_tip_check = n_tips;
  for (i=0; i<n_nodes; ++i)
  {
    pll_unode_t * t = (pll_unode_t *) malloc(sizeof(pll_unode_t));
    if (!binary_node_apply (bin_file, t, 0, bin_fread))
    {
      assert(pll_errno);
      return NULL;
    }
    if (t->next)
    {
      /* build inner node and connect */
      pll_unode_t *t_l, *t_r, *t_cl, *t_cr;
      t_l = (pll_unode_t *) malloc(sizeof(pll_unode_t));
      t_r = (pll_unode_t *) malloc(sizeof(pll_unode_t));
      retval = 1;
      retval &= binary_node_apply (bin_file, t_l, 0, bin_fread);
      retval &= binary_node_apply (bin_file, t_r, 0, bin_fread);
      if (t->label)
      {
        free(t_l->label);
        free(t_r->label);
        t_l->label = t_r->label = t->label;
      }
      if (!retval)
      {
        assert(pll_errno);
        return NULL;
      }
      t->next = t_l; t_l->next = t_r; t_r->next = t;

      /* pop */
      t_cr = tree_stack[--tree_stack_top];
      t_r->back = t_cr; t_cr->back = t_r;
      t_cl = tree_stack[--tree_stack_top];
      t_l->back = t_cl; t_cl->back = t_l;
    }
    else
      --n_tip_check;

    /* push */
    tree_stack[tree_stack_top++] = t;
  }

  /* root vertices must be in the stack */
  assert (tree_stack_top == 2);
  assert (!n_tip_check);

  tree = tree_stack[--tree_stack_top];
  tree->back = tree_stack[--tree_stack_top];
  tree->back->back = tree;

  assert(tree->pmatrix_index == tree->back->pmatrix_index);

  free(tree_stack);

  return tree;
}

/**
 *  Save a custom block to the binary file
 *
 *  @param[in] bin_file binary file
 *  @param[in] block_id id of the block for random access, or local id
 *  @param[in] data the data to store in the binary file
 *  @param[in] size the size of the data
 *  @param[in] attributes the dumped attributes
 *
 *  @return PLL_SUCCESS if the data was correctly saved
 *          PLL_FAILURE otherwise (check pll_errmsg for details)
 */
PLL_EXPORT int pllmod_binary_custom_dump(FILE * bin_file,
                                        int block_id,
                                        void * data,
                                        size_t size,
                                        unsigned int attributes)
{
  int retval;
  pll_block_header_t block_header;
  memset(&block_header, 0, sizeof(pll_block_header_t));

  /* dump header */
  block_header.block_id   = block_id;
  block_header.type       = PLLMOD_BIN_BLOCK_CUSTOM;
  block_header.attributes = attributes;
  block_header.block_len  = size;
  block_header.alignment  = 0;

  /* update main header */
  if(!binary_update_header(bin_file, &block_header))
  {
    return PLL_FAILURE;
  }

  /* dump block header */
  if (!binary_block_header_apply(bin_file, &block_header, &bin_fwrite))
    return PLL_FAILURE;

  /* dump data */
  retval = bin_fwrite(data, size, 1, bin_file);

  return retval;
}

PLL_EXPORT pll_block_map_t * pllmod_binary_get_map(FILE * bin_file,
                                                   unsigned int * n_blocks)
{
  pll_binary_header_t bin_header;
  pll_block_map_t * map;

  /* get header */
  fseek(bin_file, 0, SEEK_SET);

  if (!bin_fread(&bin_header, sizeof(pll_binary_header_t), 1, bin_file))
  {
    return NULL;
  }

  /* get map */
  map = (pll_block_map_t *) malloc (bin_header.n_blocks *
                                    sizeof(pll_block_map_t));
  if (!map) return NULL;

  if (!bin_fread(map, sizeof(pll_block_map_t), bin_header.n_blocks, bin_file))
  {
    free(map);
    return NULL;
  }

  *n_blocks = bin_header.n_blocks;
  return map;
}

/**
 *  Load a custom block from the binary file
 *
 *  @param[in] bin_file binary file
 *  @param[in] block_id id of the block for random access
 *  @param[out] size size of the block
 *  @param[out] type type of the block
 *  @param[out] attributes the block attributes
 *  @param offset offset to the data block, if known
 *                0, if access is sequential
 *                PLLMOD_BIN_ACCESS_SEEK, for searching in the file header
 *
 *  @return pointer to the loaded data
 */
PLL_EXPORT void * pllmod_binary_custom_load(FILE * bin_file,
                                           int block_id,
                                           size_t * size,
                                           unsigned int * type,
                                           unsigned int * attributes,
                                           long int offset)
{
  pll_block_header_t block_header;
  unsigned int alignment;
  void * data;

  assert (offset >= 0 || offset == PLLMOD_BIN_ACCESS_SEEK);

  if (offset != 0)
  {
    if (offset == PLLMOD_BIN_ACCESS_SEEK)
    {
      /* find offset */
      offset = binary_get_offset(bin_file, block_id);
      if (offset == PLLMOD_BIN_INVALID_OFFSET)
      {
        pllmod_set_error(PLLMOD_BIN_ERROR_BINARY_IO,
                         "Cannot retrieve offset for block %d", block_id);
        return NULL;
      }
    }

    /* apply offset */
    fseek(bin_file, offset, SEEK_SET);
  }

  /* read header */
  if (!binary_block_header_apply(bin_file, &block_header, &bin_fread))
    return NULL;
  *type       = block_header.type;
  *size       = block_header.block_len;
  alignment   = block_header.alignment;
  *attributes = block_header.attributes;

  /* read data */
  if (*attributes & PLLMOD_BIN_ATTRIB_ALIGNED)
  {
    unsigned int cur_alignment = get_current_alignment(*attributes);

    /* unimplemented so far */
    assert(cur_alignment == alignment);

    data = pll_aligned_alloc(*size, alignment);
  }
  else
    data = malloc(*size);

  if (!data)
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate space for storing data.");
    return PLL_FAILURE;
  }

  if (!bin_fread (data, *size, 1, bin_file))
  {
    pllmod_set_error(PLLMOD_BIN_ERROR_BINARY_IO, "Error reading data.");
    free(data);
    return PLL_FAILURE;
  }

  return data;
}

/* static functions */

static int cb_full_traversal(pll_unode_t * node)
{
  PLLMOD_UNUSED(node);
  return 1;
}

/**
 * Notes:
 *     1. Memory alignment could be different when saving and loading the binary
 *        file. Data will be saved without modifying the alignment, but we need
 *        to check if it is the correct one when loading.
 *     2. Binary file should be created using pllmod_binary_create. This will place
 *        the header at the beginning of the file. When a block is added, the
 *        header is updated.
 *     3. Random access binary files require some space allocated at the
 *        beginning for the hashtable.
 *
 *     4.1. Each save operation:
 *          a) Dump block header and block
 *          b) Update main header
 *          c) Update map (if necessary)
 *          d) Return file pointer to EOF
 *     4.2. Each load operation:
 *          a) If random, check map and apply offset
 *          b) Check and validate block header
 *          c) Load block
 *          d) Apply operations (e.g., new memory alignment)
  */

static unsigned int get_current_alignment( unsigned int attributes )
{
  unsigned int alignment = PLL_ALIGNMENT_CPU;
#ifdef HAVE_SSE
  if (attributes & PLL_ATTRIB_ARCH_SSE)
  alignment = PLL_ALIGNMENT_SSE;
#endif
#ifdef HAVE_AVX
  if (attributes & PLL_ATTRIB_ARCH_AVX)
  alignment = PLL_ALIGNMENT_AVX;
#endif
  return alignment;
}
