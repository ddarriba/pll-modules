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

/* bitvector */

bitv_hashtable_t *hash_init(unsigned int n,
                            unsigned int bit_count);

void hash_destroy_entry(bitv_hash_entry_t *e);

void hash_destroy(bitv_hashtable_t *h);

bitv_hash_entry_t *entry_init(double support);

hash_key_t hash_get_key(pll_split_t s, int len);

int hash_update(pll_split_t bit_vector,
                 bitv_hashtable_t *h,
                 hash_key_t key,
                 double support,
                 unsigned int position);

void hash_insert(pll_split_t bit_vector,
                 bitv_hashtable_t *h,
                 unsigned int bip_number,
                 hash_key_t key,
                 double support,
                 unsigned int position);

void hash_remove(bitv_hashtable_t *h,
                 bitv_hash_entry_t ** prev_ptr,
                 bitv_hash_entry_t * e);

void hash_print(bitv_hashtable_t *h);

/* bitvector utilities */

void bitv_normalize(pll_split_t bitv, unsigned int bit_count);

int bitv_is_normalized(const pll_split_t bitv);

unsigned int bitv_length(unsigned int bit_count);

/* string */

string_hashtable_t *string_hash_init(unsigned int n, unsigned int max_labels);

void string_hash_destroy(string_hashtable_t *h);

hash_key_t string_hash_get_key(const char * s);

int string_hash_insert(const char *s,
                       string_hashtable_t *h,
                       int node_number);

int string_hash_lookup(char *s, string_hashtable_t *h);

#endif
