/*
 Copyright (C) 2016 Diego Darriba, Alexandros Stamatakis

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

#include "tree_hashtable.h"
#include "../pllmod_common.h"

bitv_hashtable_t *hash_init(unsigned int n,
                            unsigned int bit_count)
{
  /* init with powers of two */
  static const unsigned int init_table[] = {
                64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384,
					      32768, 65536, 131072, 262144, 524288, 1048576, 2097152,
					      4194304, 8388608, 16777216, 33554432, 67108864, 134217728,
					      268435456, 536870912, 1073741824, 2147483648U};

  bitv_hashtable_t *h = (bitv_hashtable_t*) malloc(sizeof(bitv_hashtable_t));
  if (!h)
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for hashtable\n");
    return NULL;
  }

  unsigned int
    table_size,
    i,
    init_table_size = sizeof(init_table)/sizeof(init_table[0]),
    maxSize = (unsigned int)-1;

  assert(n <= maxSize);

  i = 0;

  while(init_table[i] < n && i < init_table_size)
    ++i;

  assert(i < init_table_size);

  table_size = init_table[i];

  /* printf("Hash table init with size %u\n", table_size); */

  h->table = (bitv_hash_entry_t**)calloc(table_size, sizeof(bitv_hash_entry_t*));
  h->table_size = table_size;
  h->entry_count = 0;
  h->bit_count = bit_count;
  h->bitv_len = bitv_length(bit_count);

  if (!h->table)
  {
    free(h);
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for hashtable entries\n");
    return NULL;
  }

  return h;
}

void hash_destroy_entry(bitv_hash_entry_t *e)
{
  if(e->bit_vector)
    free(e->bit_vector);

  if(e->tree_vector)
    free(e->tree_vector);

  free(e);
}

void hash_destroy(bitv_hashtable_t *h)
{
  unsigned int
    i,
    entry_count = 0;


  for(i = 0; i < h->table_size; ++i)
  {
    if(h->table[i] != NULL)
  	{
  	  bitv_hash_entry_t *e = h->table[i];
  	  bitv_hash_entry_t *previous;

  	  do
	    {
	      previous = e;
	      e = e->next;

        hash_destroy_entry(previous);
	      ++entry_count;
	    }
  	  while(e != NULL);
  	}
  }

  assert(entry_count == h->entry_count);

  free(h->table);
  free(h);
}

bitv_hash_entry_t *entry_init(double support)
{
  bitv_hash_entry_t *e = (bitv_hash_entry_t*)malloc(sizeof(bitv_hash_entry_t));

  e->bit_vector     = (pll_split_t)NULL;
  e->tree_vector    = (unsigned int*)NULL;
  e->support        = support;
  e->bip_number     = 0;
  e->next       = (bitv_hash_entry_t*)NULL;

  return e;
}

hash_key_t hash_get_key(pll_split_t s, int len)
{
  hash_key_t h = 0;
  int i;

  for(i = 0; i < len; ++i)
  {
    h += s[i];
    h += ( h << 10 );
    h ^= ( h >> 6 );
  }

  h += ( h << 3 );
  h ^= ( h >> 11 );
  h += ( h << 15 );

  return h;
}

/* this function only increments support for existing splits,
 * but never adds new splits to the hashtable */
int hash_update(pll_split_t bit_vector,
                bitv_hashtable_t *h,
                hash_key_t key,
                double support,
                unsigned int position)
{
  if (key == HASH_KEY_UNDEF)
  {
      key = hash_get_key(bit_vector, (int)h->bitv_len);
      position = key % h->table_size;
  }

  if(h->table[position] != NULL)
  {
    bitv_hash_entry_t *e = h->table[position];
    do
    {
      unsigned int i = 0;

      /* check for identity of bipartitions */

      if (e->key == key)
        for(i = 0; i < h->bitv_len; ++i)
          if(bit_vector[i] != e->bit_vector[i])
            break;

      if(i == h->bitv_len)
      {
        e->support = e->support + support;
        return PLL_SUCCESS;
      }

      /* otherwise keep searching */
      e = e->next;
    }
    while(e != (bitv_hash_entry_t*)NULL);
  }

  return PLL_FAILURE;
}

void hash_insert(pll_split_t bit_vector,
                 bitv_hashtable_t *h,
                 unsigned int bip_number,
                 hash_key_t key,
                 double support,
                 unsigned int position)
{
  bitv_hash_entry_t *e;

  if (key == HASH_KEY_UNDEF)
  {
      key = hash_get_key(bit_vector, (int)h->bitv_len);
      position = key % h->table_size;
  }

  if(h->table[position] != NULL)
    {
      /* search for this split in hashtable, and increment its support if found */
      if (hash_update(bit_vector, h, key, support, position))
        return;

      /* add new split to the hashtable */
      e = entry_init(support);
      e->key = key;
      e->bip_number = bip_number;

      e->bit_vector = (pll_split_t) calloc(h->bitv_len, sizeof(pll_split_base_t));
      memcpy(e->bit_vector, bit_vector, sizeof(pll_split_base_t) * h->bitv_len);

      e->next = h->table[position];
      h->table[position] = e;
    }
  else
  {
    e = entry_init(support);
    e->key = key;
    e->bip_number = bip_number;
    e->bit_vector = (pll_split_t) calloc(h->bitv_len, sizeof(pll_split_base_t));
    memcpy(e->bit_vector, bit_vector, sizeof(pll_split_base_t) * h->bitv_len);

    h->table[position] = e;
  }

  h->entry_count =  h->entry_count + 1;
}

void hash_remove(bitv_hashtable_t *h,
                 bitv_hash_entry_t ** prev_ptr,
                 bitv_hash_entry_t * e)
{
  *prev_ptr = e->next;
  hash_destroy_entry(e);
  --h->entry_count;
  assert(h->entry_count >= 0);
}

void hash_print(bitv_hashtable_t *h)
{
  unsigned int i;
  for (i=0; i<h->table_size; ++i)
  {
    bitv_hash_entry_t * e =  h->table[i];
    while (e != NULL)
    {
      pllmod_utree_split_show(e->bit_vector, h->bit_count);
      printf(" %f\n", e->support);
      e = e->next;
    }
  }
}

void bitv_normalize(pll_split_t bitv, unsigned int bit_count)
{
  unsigned int split_size  = sizeof(pll_split_base_t) * 8;
  unsigned int split_offset = bit_count % split_size;
  unsigned int split_len    = bitv_length(bit_count);
  unsigned int i;

  int normalized = bitv_is_normalized(bitv);

  if (!normalized)
  {
    for (i=0; i<split_len; ++i)
    {
      bitv[i] = ~bitv[i];
    }

    if (split_offset)
    {
      unsigned int mask = (1<<split_offset) - 1;
      bitv[split_len - 1] &= mask;
    }
  }
}

int bitv_is_normalized(const pll_split_t bitv)
{
  return bitv[0]&1;
}

unsigned int bitv_length(unsigned int bit_count)
{
  unsigned int split_size = sizeof(pll_split_base_t) * 8;
  unsigned int split_offset = bit_count % split_size;

  return bit_count / split_size + (split_offset>0);
}

/* string */

string_hashtable_t *string_hash_init(unsigned int n, unsigned int max_labels)
{
  /* init with primes */
  static const hash_key_t init_table[] = {
                                  53, 97, 193, 389, 769, 1543, 3079, 6151,
                                  12289, 24593, 49157, 98317, 196613, 393241,
                                  786433, 1572869, 3145739, 6291469, 12582917,
                                  25165843, 50331653, 100663319, 201326611,
                                  402653189, 805306457, 1610612741};

  string_hashtable_t *h = (string_hashtable_t*)malloc(sizeof(string_hashtable_t));
  assert(h);

  unsigned int
    table_size,
    i,
    prime_table_length = sizeof(init_table)/sizeof(hash_key_t),
    max_size = (hash_key_t)-1;

  assert(n <= max_size);

  i = 0;

  while(init_table[i] < n && i < prime_table_length)
    ++i;

  assert(i < prime_table_length);

  table_size = init_table[i];

  h->table = (string_hash_entry_t**)calloc(table_size,
                                           sizeof(string_hash_entry_t*));
  h->labels = (char **) malloc(max_labels * sizeof(char *));
  h->table_size = table_size;
  h->entry_count = 0;

  return h;
}

void string_hash_destroy(string_hashtable_t *h)
{
  unsigned int entry_count = 0, i;

  for(i = 0; i < h->table_size; ++i)
  {
    if(h->table[i] != NULL)
  	{
  	  string_hash_entry_t *e = h->table[i];
  	  string_hash_entry_t *previous;

  	  do
	    {
	      previous = e;
	      e = e->next;

	      if(previous->word)
		      free(previous->word);

	      free(previous);
	      ++entry_count;
	    }
  	  while(e != NULL);
  	}
  }

  assert(entry_count == h->entry_count);

  free(h->labels);
  free(h->table);
  free(h);
}

hash_key_t string_hash_get_key(const char * s)
{
  hash_key_t h = 0;

  for(; *s; ++s)
    h = 31 * h + (unsigned int)*s;

  return h;
}

int string_hash_insert(const char *s,
                string_hashtable_t *h,
                int node_number)
{
  hash_key_t position = string_hash_get_key(s) % h->table_size;
  string_hash_entry_t *p = h->table[position];

  for(; p!= NULL; p = p->next)
  {
    if(strcmp(s, p->word) == 0)
	    return PLL_FAILURE;
  }

  p = (string_hash_entry_t *)malloc(sizeof(string_hash_entry_t));
  assert(p);

  h->labels[node_number] = (char *) malloc((strlen(s) + 1) * sizeof(char));
  p->node_number = node_number;
  p->word = h->labels[node_number];
  strcpy(p->word, s);

  p->next =  h->table[position];
  h->table[position] = p;
  ++h->entry_count;

  return PLL_SUCCESS;
}

int string_hash_lookup(char *s, string_hashtable_t *h)
{
  hash_key_t position = string_hash_get_key(s) % h->table_size;
  string_hash_entry_t *p = h->table[position];

  for(; p!= NULL; p = p->next)
  {
    if(strcmp(s, p->word) == 0)
      return p->node_number;
  }

  return -1;
}
