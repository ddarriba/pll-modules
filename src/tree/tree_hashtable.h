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

#ifndef TREE_HASHTABLE_H_
#define TREE_HASHTABLE_H_

#include "pll_tree.h"

typedef unsigned int hash_key_t;

typedef struct bitv_hash_entry
{
  hash_key_t key;
  pll_split_t bit_vector;
  unsigned int *tree_vector;
  unsigned int tip_count;
  int support;
  unsigned int bip_number;

  struct bitv_hash_entry *next;
} bitv_hash_entry_t;

typedef struct
{
  unsigned int table_size;
  bitv_hash_entry_t **table;
  unsigned int entry_count;
} bitv_hashtable_t;

typedef struct string_hash_entry
{
  hash_key_t key;
  int node_number;
  char * word;
  struct string_hash_entry *next;
} string_hash_entry_t;

typedef struct
{
  char **labels;
  unsigned int table_size;
  string_hash_entry_t **table;
  unsigned int entry_count;
} string_hashtable_t;

/* bitvector */

bitv_hashtable_t *hash_init(unsigned int n);

void hash_destroy(bitv_hashtable_t *h);

bitv_hash_entry_t *entry_init(void);

hash_key_t hash_get_key(pll_split_t s, int len);

void hash_insert(pll_split_t bit_vector,
                 bitv_hashtable_t *h,
                 unsigned int vector_length,
                 unsigned int bip_number,
                 hash_key_t key,
                 unsigned int position);

/* string */

string_hashtable_t *string_hash_init(unsigned int n, unsigned int max_labels);

void string_hash_destroy(string_hashtable_t *h);

hash_key_t string_hash_get_key(const char * s);

int string_hash_insert(const char *s,
                       string_hashtable_t *h,
                       int node_number);

int string_hash_lookup(char *s, string_hashtable_t *h);

#endif
