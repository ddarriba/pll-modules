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
#include "binary_io_operations.h"

/**
 * Notes:
 *     1. Memory alignment could be different when saving and loading the binary
 *        file. Data will be saved without modifying the alignment, but we need
 *        to check if it is the correct one when loading.
 *     2. Binary file should be created using pll_binary_create. This will place
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

static unsigned long clv_size(pll_partition_t * partition);
static unsigned long partition_size(pll_partition_t * partition);
static unsigned long weightv_size(pll_partition_t * partition);
static unsigned long scaler_size(pll_partition_t * partition);

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

static void file_io_error (FILE * bin_file, long int setp, const char * msg)
{
  assert(setp >= PLL_BINARY_INVALID_OFFSET);

  /* if offset is valid, we apply it */
  if (setp != PLL_BINARY_INVALID_OFFSET)
    fseek(bin_file, setp, SEEK_SET);

  /* update error data */
  snprintf(pll_errmsg, 200, "Binary file I/O error: %s", msg);
  pll_errno = PLL_ERROR_LOADSTORE;
}

static int binary_update_header(FILE * bin_file,
                               pll_block_header_t * header)
{
  unsigned int next_block;
  pll_block_map_t next_map;

  long int cur_position = ftell(bin_file);

  pll_binary_header_t bin_header;

  /* update header */
  fseek(bin_file, 0, SEEK_SET);
  if (!bin_fread(&bin_header, sizeof(pll_binary_header_t), 1, bin_file))
  {
    file_io_error(bin_file, cur_position, "update binary header [r]");
    return PLL_FAILURE;
  }
  fseek(bin_file, 0, SEEK_SET);
  next_block = bin_header.n_blocks;
  ++bin_header.n_blocks;

  if (!bin_fwrite(&bin_header, sizeof(pll_binary_header_t), 1, bin_file))
  {
    file_io_error(bin_file, cur_position, "update binary header [w]");
    return PLL_FAILURE;
  }

  if (header && (header->attributes & PLL_BINARY_ATTRIB_UPDATE_MAP))
  {
    /* update map */
    assert(next_block < bin_header.max_blocks);
    fseek(bin_file, next_block * sizeof(pll_block_map_t), SEEK_CUR);

    next_map.block_id     = header->block_id;
    next_map.block_offset = cur_position;
    if (!bin_fwrite(&next_map, sizeof(pll_block_map_t), 1, bin_file))
    {
      file_io_error(bin_file, cur_position, "update binary map");
      return PLL_FAILURE;
    }
  }

  /* move back to original position */
  fseek(bin_file, cur_position, SEEK_SET);
  return PLL_SUCCESS;
}

PLL_EXPORT FILE * pll_binary_create(const char * filename,
                                    pll_binary_header_t * header,
                                    unsigned int access_type,
                                    unsigned int n_blocks)
{
  FILE * file = NULL;

  memset(header, 0, sizeof(pll_binary_header_t));
  header->access_type = access_type;
  header->max_blocks = n_blocks;
  header->map_offset = (access_type == PLL_BINARY_ACCESS_RANDOM)?
      n_blocks * sizeof(pll_block_map_t):0;
  header->n_blocks = 0;

  if (access_type == PLL_BINARY_ACCESS_RANDOM && n_blocks <= 0)
  {
    snprintf(pll_errmsg, 200,
             "Number of blocks for random access must be greater than 0");
    pll_errno = PLL_ERROR_INVALID_SIZE;
    return NULL;
  }

  file = fopen(filename, "w+b");

  if (!file)
  {
    snprintf(pll_errmsg, 200, "Cannot open file for writing");
    pll_errno = PLL_ERROR_FILE_OPEN;
    return NULL;
  }

  if (!bin_fwrite(header, sizeof(pll_binary_header_t), 1, file))
  {
    snprintf(pll_errmsg, 200, "Error writing header to file");
    pll_errno = PLL_ERROR_BINARY_IO;
    fclose(file);
    return NULL;
  }

  if(fseek(file, header->map_offset, SEEK_CUR))
  {
    snprintf(pll_errmsg, 200, "Error seeking through file during creation");
    pll_errno = PLL_ERROR_BINARY_IO;
    fclose(file);
    return NULL;
  }

  return file;
}

PLL_EXPORT FILE * pll_binary_open(const char * filename,
                                  pll_binary_header_t * header)
{
  FILE * file;

  file = fopen(filename, "r");

  if (!file)
  {
    snprintf(pll_errmsg, 200, "Cannot open file for reading");
    pll_errno = PLL_ERROR_FILE_OPEN;
    return NULL;
  }

  if (!bin_fread(header, sizeof(pll_binary_header_t), 1, file))
  {
    snprintf(pll_errmsg, 200, "Error reading header from file");
    pll_errno = PLL_ERROR_BINARY_IO;
    fclose(file);
    return NULL;
  }

  fseek(file, header->map_offset, SEEK_CUR);

  return file;
}

PLL_EXPORT int pll_binary_close(FILE * bin_file)
{
  return fclose(bin_file);
}

PLL_EXPORT int pll_binary_partition_dump(FILE * bin_file,
                                         int block_id,
                                         pll_partition_t * partition,
                                         unsigned int attributes)
{
  pll_block_header_t block_header;
  unsigned long partition_len = partition_size(partition),
                clv_len = 0,
                wgt_len = 0;

  /* fill block header */
  block_header.block_id   = block_id;
  block_header.type       = PLL_BINARY_BLOCK_PARTITION;
  block_header.attributes = attributes;
  block_header.block_len  = partition_len;
  block_header.alignment  = 0;

  if (attributes & PLL_BINARY_ATTRIB_PARTITION_DUMP_CLV)
  {
    clv_len = (partition->clv_buffers + partition->tips)  * clv_size(partition) +
               partition->scale_buffers * scaler_size(partition);
    block_header.block_len += clv_len;
  }
  if (attributes & PLL_BINARY_ATTRIB_PARTITION_DUMP_WGT)
    {
      wgt_len = weightv_size(partition);
      block_header.block_len += wgt_len;
    }

  /* update main header */
  binary_update_header(bin_file, &block_header);

  /* dump header */
  bin_fwrite(&block_header, sizeof(pll_block_header_t) ,1, bin_file);

  /* dump data */
  binary_apply_to_partition(bin_file, partition, attributes, &bin_fwrite);

  return PLL_SUCCESS;
}

PLL_EXPORT pll_partition_t * pll_binary_partition_load(FILE * bin_file,
                                                       int block_id,
                                                       pll_partition_t * partition,
                                                       unsigned int * attributes,
                                                       const unsigned int * map,
                                                       long int offset)
{
  pll_block_header_t block_header;
  pll_partition_t * local_partition;
  assert(offset >= 0 || offset == PLL_BINARY_ACCESS_SEEK);

  if (offset != 0)
  {
    if (offset == PLL_BINARY_ACCESS_SEEK)
    {
      /* find offset */
      offset = binary_get_offset (bin_file, block_id);
      if (offset == PLL_BINARY_INVALID_OFFSET)
        return NULL;
    }

    /* apply offset */
    fseek (bin_file, offset, SEEK_SET);
  }

  bin_fread (&block_header, sizeof(pll_block_header_t), 1, bin_file);

  if (block_header.type != PLL_BINARY_BLOCK_PARTITION)
  {
    pll_errno = PLL_ERROR_BLOCK_MISMATCH;
    snprintf(pll_errmsg, 200, "Block type is %d and should be %d",
             block_header.type, PLL_BINARY_BLOCK_PARTITION);
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
    binary_apply_to_partition_desc (bin_file, &aux_partition, *attributes, &bin_fread);
    local_partition = pll_partition_create(
        aux_partition.tips,
        aux_partition.clv_buffers,
        aux_partition.states,
        aux_partition.sites,
        aux_partition.rate_matrices,
        aux_partition.prob_matrices,
        aux_partition.rate_cats,
        aux_partition.scale_buffers,
        map,
        aux_partition.attributes);

    if (!local_partition)
      return NULL;
  }

  binary_apply_to_partition_body (bin_file, local_partition, *attributes, &bin_fread);

  return local_partition;
}


PLL_EXPORT int pll_binary_clv_dump(FILE * bin_file,
                                   int block_id,
                                   pll_partition_t * partition,
                                   unsigned int clv_index,
                                   unsigned int attributes)
{
  int retval;
  pll_block_header_t block_header;

  size_t clv_size = partition->sites * partition->states_padded *
                      partition->rate_cats;

  /* fill block header */
  block_header.block_id   = block_id;
  block_header.type       = PLL_BINARY_BLOCK_CLV;
  block_header.attributes = attributes;
  block_header.block_len  = clv_size * sizeof(double);
  block_header.alignment  = 0;

  /* update main header */
  binary_update_header(bin_file, &block_header);

  /* dump block header */
  bin_fwrite(&block_header, sizeof(pll_block_header_t) ,1, bin_file);

  /* dump data */
  retval = binary_apply_to_clv (bin_file,
                         partition,
                         clv_index,
                         attributes,
                         clv_size,
                         &bin_fwrite);

  return retval;
}

PLL_EXPORT int pll_binary_clv_load(FILE * bin_file,
                                   int block_id,
                                   pll_partition_t * partition,
                                   unsigned int clv_index,
                                   unsigned int * attributes,
                                   long int offset)
{
  int retval;
  pll_block_header_t block_header;

  assert (partition);
  assert(offset >= 0 || offset == PLL_BINARY_ACCESS_SEEK);

  if (offset != 0)
  {
    if (offset == PLL_BINARY_ACCESS_SEEK)
    {
      /* find offset */
      offset = binary_get_offset (bin_file, block_id);
      if (offset == PLL_BINARY_INVALID_OFFSET)
        return PLL_FAILURE;
    }

    /* apply offset */
    fseek (bin_file, offset, SEEK_SET);
  }

  /* read and validate header */
  bin_fread (&block_header, sizeof(pll_block_header_t), 1, bin_file);
  if (block_header.type != PLL_BINARY_BLOCK_CLV)
  {
    pll_errno = PLL_ERROR_BLOCK_MISMATCH;
    snprintf(pll_errmsg, 200, "Block type is %d and should be %d",
             block_header.type, PLL_BINARY_BLOCK_CLV);
    return PLL_FAILURE;
  }

  size_t clv_size = partition->sites * partition->states_padded *
                    partition->rate_cats;

  if (block_header.block_len != (clv_size * sizeof(double)))
  {
      pll_errno = PLL_ERROR_BLOCK_LENGTH;
      snprintf(pll_errmsg, 200, "Wrong block length");
      return PLL_FAILURE;
  }

  *attributes = block_header.attributes;

  retval = binary_apply_to_clv (bin_file,
                         partition,
                         clv_index,
                         *attributes,
                         clv_size,
                         &bin_fread);

  return retval;
}

static int cb_full_traversal(pll_utree_t * node)
{
  return 1;
}

PLL_EXPORT int pll_binary_utree_dump(FILE * bin_file,
                                     int block_id,
                                     pll_utree_t * tree,
                                     unsigned int tip_count,
                                     unsigned int attributes)
{
  pll_utree_t ** travbuffer;
  pll_block_header_t block_header;
  unsigned int i, n_nodes, n_inner, n_utrees, trav_size;

  n_inner = tip_count - 2;
  n_nodes = tip_count + n_inner;
  n_utrees = tip_count + 3 * n_inner;

  travbuffer = (pll_utree_t **)malloc(n_nodes* sizeof(pll_utree_t *));

  pll_utree_traverse(tree, cb_full_traversal, travbuffer, &trav_size);

  assert (trav_size == n_nodes);

  block_header.block_id   = block_id;
  block_header.type       = PLL_BINARY_BLOCK_TREE;
  block_header.attributes = attributes;
  block_header.block_len  = n_utrees * sizeof(pll_utree_t);
  block_header.alignment  = 0;

  /* update main header */
  if(!binary_update_header(bin_file, &block_header))
  {
    return PLL_FAILURE;
  }

  /* dump block header */
  bin_fwrite(&block_header, sizeof(pll_block_header_t) ,1, bin_file);

  /* dump data */
  for (i=0; i<trav_size; ++i)
  {
    binary_apply_to_node (bin_file, travbuffer[i], 1, bin_fwrite);
    if (travbuffer[i]->next)
    {
      binary_apply_to_node (bin_file, travbuffer[i]->next, 1, bin_fwrite);
      binary_apply_to_node (bin_file, travbuffer[i]->next->next, 1, bin_fwrite);
    }
  }

  free(travbuffer);

  return PLL_SUCCESS;
}

PLL_EXPORT pll_utree_t * pll_binary_utree_load(FILE * bin_file,
                                               int block_id,
                                               unsigned int * attributes,
                                               long int offset)
{
  unsigned int i, n_tips, n_tip_check, n_nodes;
  long n_utrees;
  pll_block_header_t block_header;
  pll_utree_t ** tree_stack;
  pll_utree_t * tree;
  unsigned int tree_stack_top;

  assert(offset >= 0 || offset == PLL_BINARY_ACCESS_SEEK);

  if (offset != 0)
  {
    if (offset == PLL_BINARY_ACCESS_SEEK)
    {
      /* find offset */
      offset = binary_get_offset (bin_file, block_id);
      if (offset == PLL_BINARY_INVALID_OFFSET)
        return PLL_FAILURE;
    }

    /* apply offset */
    fseek (bin_file, offset, SEEK_SET);
  }

  /* read and validate header */
  bin_fread (&block_header, sizeof(pll_block_header_t), 1, bin_file);
  if (block_header.type != PLL_BINARY_BLOCK_TREE)
  {
    pll_errno = PLL_ERROR_BLOCK_MISMATCH;
    snprintf (pll_errmsg, 200, "Block type is %d and should be %d",
              block_header.type, PLL_BINARY_BLOCK_TREE);
    return PLL_FAILURE;
  }

  n_utrees = block_header.block_len/sizeof(pll_utree_t);
  n_tips   = (unsigned int) ((n_utrees + 6) / 4);
  n_nodes  = 2*n_tips - 2;
  assert( n_utrees % 4 == 2 );

  *attributes = block_header.attributes;

  /* allocate stack for at most 'n_tips' nodes */
  tree_stack = (pll_utree_t **) malloc(n_tips * sizeof (pll_utree_t *));
  tree_stack_top = 0;

  /* read nodes */
  n_tip_check = n_tips;
  for (i=0; i<n_nodes; ++i)
  {
    pll_utree_t * t = (pll_utree_t *) malloc(sizeof(pll_utree_t));
    binary_apply_to_node (bin_file, t, 0, bin_fread);
    if (t->next)
    {
      /* build inner node and connect */
      pll_utree_t *t_l, *t_r, *t_cl, *t_cr;
      t_l = (pll_utree_t *) malloc(sizeof(pll_utree_t));
      t_r = (pll_utree_t *) malloc(sizeof(pll_utree_t));
      binary_apply_to_node (bin_file, t_l, 0, bin_fread);
      binary_apply_to_node (bin_file, t_r, 0, bin_fread);
      if (t->label)
      {
        free(t_l->label);
        free(t_r->label);
        t_l->label = t_r->label = t->label;
      }
      t->next = t_l; t_l->next = t_r; t_r->next = t;

      assert(t->clv_index == t_l->clv_index);
      assert(t->clv_index == t_r->clv_index);
      assert(t->scaler_index == t_l->scaler_index);
      assert(t->scaler_index == t_r->scaler_index);

      /* pop */
      t_cl = tree_stack[--tree_stack_top];
      t_l->back = t_cl; t_cl->back = t_l;
      assert(t_l->pmatrix_index == t_cl->back->pmatrix_index);
      t_cr = tree_stack[--tree_stack_top];
      t_r->back = t_cr; t_cr->back = t_r;
      assert(t_r->pmatrix_index == t_cr->back->pmatrix_index);
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

PLL_EXPORT int pll_binary_custom_dump(FILE * bin_file,
                                      int block_id,
                                      void * data,
                                      size_t size,
                                      unsigned int attributes)
{
  int retval;
  pll_block_header_t block_header;
  memset(&block_header, 0, sizeof(pll_block_header_t));

  /* update main header */
  binary_update_header(bin_file, &block_header);

  /* dump header */
  block_header.block_id   = block_id;
  block_header.type       = PLL_BINARY_BLOCK_CLV;
  block_header.attributes = attributes;
  block_header.block_len  = size;
  block_header.alignment  = 0;
  bin_fwrite(&block_header, sizeof(pll_block_header_t) ,1, bin_file);

  /* dump data */
  retval = bin_fwrite(data, size, 1, bin_file);

  return retval;
}

PLL_EXPORT pll_block_map_t * pll_binary_get_map(FILE * bin_file,
                                                unsigned int * n_blocks)
{
  pll_binary_header_t bin_header;
  pll_block_map_t * map;

  /* get header */
  fseek(bin_file, 0, SEEK_SET);

  if (!bin_fread(&bin_header, sizeof(pll_binary_header_t), 1, bin_file))
  {
    file_io_error(bin_file, PLL_BINARY_INVALID_OFFSET, "read binary header");
    return NULL;
  }

  /* get map */
  map = (pll_block_map_t *) malloc (bin_header.n_blocks * sizeof(pll_block_map_t));
  if (!map) return NULL;

  if (!bin_fread(map, sizeof(pll_block_map_t), bin_header.n_blocks, bin_file))
  {
    free(map);
    file_io_error(bin_file, PLL_BINARY_INVALID_OFFSET, "read binary header");
    return NULL;
  }

  *n_blocks = bin_header.n_blocks;
  return map;
}

PLL_EXPORT void * pll_binary_custom_load(FILE * bin_file,
                                         int block_id,
                                         size_t * size,
                                         unsigned int * type,
                                         unsigned int * attributes,
                                         long int offset)
{
  pll_block_header_t block_header;
  unsigned int alignment;
  void * data;

  assert (offset >= 0 || offset == PLL_BINARY_ACCESS_SEEK);

  if (offset != 0)
  {
    if (offset == PLL_BINARY_ACCESS_SEEK)
    {
      /* find offset */
      offset = binary_get_offset(bin_file, block_id);
      if (offset == PLL_BINARY_INVALID_OFFSET)
        return NULL;
    }

    /* apply offset */
    fseek(bin_file, offset, SEEK_SET);
  }

  /* read header */
  bin_fread (&block_header, sizeof(pll_block_header_t), 1, bin_file);
  *type       = block_header.type;
  *size       = block_header.block_len;
  alignment   = block_header.alignment;
  *attributes = block_header.attributes;

  /* read data */
  if (*attributes & PLL_BINARY_ATTRIB_ALIGNED)
  {
    unsigned int cur_alignment = get_current_alignment(*attributes);

    /* unimplemented so far */
    assert(cur_alignment == alignment);

    data = pll_aligned_alloc(*size,alignment);
  }
  else
    data = malloc(*size);

  if (!data)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf (pll_errmsg, 200,
              "Cannot allocate space for storing data.");
    return PLL_FAILURE;
  }

  if (!bin_fread (data, *size, 1, bin_file))
  {
    pll_errno = PLL_ERROR_BINARY_IO;
    snprintf (pll_errmsg, 200, "Error reading data.");
    free(data);
    return PLL_FAILURE;
  }

  return data;
}

static unsigned long clv_size(pll_partition_t * partition)
{
  unsigned long size = 0;
  unsigned int n_clvs, n_tipchars;

  /* Question: What will happen when tipchars are set? partition->clv_buffers
   * would count the tips as well even though they are not used?
   */
  n_tipchars = partition->attributes & PLL_ATTRIB_PATTERN_TIP?partition->tips:0;
  n_clvs     = partition->clv_buffers;

  size = sizeof(double) *
         partition->sites *
         partition->states_padded *
         partition->rate_cats *
         n_clvs;

  if (n_tipchars)
  {
    size += sizeof(char) * partition->sites * n_tipchars;
    size += sizeof(char) * PLL_ASCII_SIZE;
  }

  return size;
}

static unsigned long scaler_size(pll_partition_t * partition)
{
  unsigned long size = 0;
  /*TODO: when per-category scalers are available, we need to multiply
   * the size by the number of categories.
   */
  size = sizeof(unsigned int)    *
         partition->sites;
  return size;
}

static unsigned long weightv_size(pll_partition_t * partition)
{
  unsigned long size = 0;
  size = sizeof(unsigned int)    *
         partition->sites;
  return size;
}

static unsigned long partition_size(pll_partition_t * partition)
{
  unsigned long size         = 0;
  unsigned int rate_cats     = partition->rate_cats;
  unsigned int states        = partition->states;
  unsigned int states_padded = partition->states_padded;
  unsigned int n_subst_rates = states * (states-1) / 2;
  unsigned int prob_matrices = partition->prob_matrices;
  unsigned int rate_matrices = partition->rate_matrices;

//  tips clv_buffers states sites rate_matrices prob_matrices rate_cats;
//  scale_buffers attributes states_padded maxstates
//  log2_maxstates log2_states log2_rates
  size += sizeof(unsigned int) * 15;

  /* eigen_decomp_valid */
  size += sizeof(int) * rate_matrices;
  /* eigenvecs + inv_eigenvecs */
  size += 2 * (sizeof(double) * rate_matrices *
               states_padded * states_padded);
  /* eigenvals */
  size += sizeof(double) * rate_matrices * states_padded;

  /* pmatrix */
  size += sizeof(double) * prob_matrices * states_padded * states_padded * rate_cats;

  /* subst_params */
  size += sizeof(double) * rate_matrices * n_subst_rates;

  /* frequencies */
  size += sizeof(double) * rate_matrices * states_padded;

  /* rate_cats + rate_weights */
  size += 2 * sizeof(double) * rate_cats;

  /* prop_invar */
  size += 2 * sizeof(double) * rate_matrices;

  return size;
}
