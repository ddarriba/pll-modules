/*
 * binary_io_operations.h
 *
 *  Created on: Mar 23, 2016
 *      Author: diego
 */
#ifndef BINARY_IO_OPERATIONS_H_
#define BINARY_IO_OPERATIONS_H_

#include "pll_binary.h"

int bin_fread(void * data, size_t size, size_t count, FILE * file);
int bin_fwrite(void * data, size_t size, size_t count, FILE * file);

long int binary_get_offset(FILE *bin_file, int block_id);

int binary_apply_to_partition(FILE * bin_file,
                       pll_partition_t * partition,
                       unsigned int attributes,
                       int (*bin_func)(void *, size_t, size_t, FILE *));

int binary_apply_to_partition_body (FILE * bin_file,
                             pll_partition_t * partition,
                             unsigned int attributes,
                             int (*bin_func)(void *, size_t, size_t, FILE *));

int binary_apply_to_partition_desc (FILE * bin_file,
                             pll_partition_t * partition,
                             unsigned int attributes,
                             int (*bin_func)(void *, size_t, size_t, FILE *));

int binary_apply_to_clv (FILE * bin_file,
                  pll_partition_t * partition,
                  unsigned int clv_index,
                  unsigned int attributes,
                  size_t clv_size,
                  int (*bin_func)(void *, size_t, size_t, FILE *));

int binary_apply_to_node (FILE * bin_file,
                          pll_utree_t * node,
                          int write,
                          int (*bin_func)(void *, size_t, size_t, FILE *));

#endif /* BINARY_IO_OPERATIONS_H_ */
