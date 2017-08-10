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
#ifndef BINARY_IO_OPERATIONS_H_
#define BINARY_IO_OPERATIONS_H_

#include "pll_binary.h"

int bin_fread(void * data, size_t size, size_t count, FILE * file);

int bin_fwrite(void * data, size_t size, size_t count, FILE * file);

int binary_block_header_apply(FILE * bin_file,
                              pll_block_header_t * block_header,
                              int (*bin_func)(void *, size_t, size_t, FILE *));

int binary_update_header(FILE * bin_file,
                         pll_block_header_t * header);

long int binary_get_offset(FILE *bin_file, int block_id);

int binary_partition_apply(FILE * bin_file,
                           pll_partition_t * partition,
                           unsigned int attributes,
                           int (*bin_func)(void *, size_t, size_t, FILE *));

int binary_partition_body_apply (FILE * bin_file,
                             pll_partition_t * partition,
                             unsigned int attributes,
                             int (*bin_func)(void *, size_t, size_t, FILE *));

int binary_partition_desc_apply (FILE * bin_file,
                             pll_partition_t * partition,
                             unsigned int attributes,
                             int (*bin_func)(void *, size_t, size_t, FILE *));

int binary_repeats_apply (FILE * bin_file,
                  pll_partition_t * partition,
                  unsigned int attributes,
                  size_t nodes,
                  int (*bin_func)(void *, size_t, size_t, FILE *));

int binary_clv_apply (FILE * bin_file,
                  pll_partition_t * partition,
                  unsigned int clv_index,
                  unsigned int attributes,
                  size_t clv_size,
                  int (*bin_func)(void *, size_t, size_t, FILE *));

int binary_node_apply (FILE * bin_file,
                       pll_unode_t * node,
                       int write,
                       int (*bin_func)(void *, size_t, size_t, FILE *));

void file_io_error (FILE * bin_file, long int setp, const char * msg);

#endif /* BINARY_IO_OPERATIONS_H_ */
