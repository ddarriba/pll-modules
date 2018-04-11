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
  * @file utree_distances.c
  *
  * @brief Operations on splits for unrooted trees
  *
  * @author Diego Darriba
  */

#include "pll_tree.h"
#include "tree_hashtable.h"

#include "../pllmod_common.h"

static int cb_get_splits(pll_unode_t * node, void *data);
static inline void merge_split(pll_split_t to,
                         const pll_split_t from,
                         unsigned int split_len);
static int _cmp_splits (const void * a, const void * b);
static int _cmp_split_node_pair (const void * a, const void * b);
static int compare_splits (pll_split_t s1,
                           pll_split_t s2,
                           unsigned int split_len);
static unsigned int get_utree_splitmap_id(pll_unode_t * node,
                                          unsigned int tip_count);
static int split_is_valid_and_normalized(const pll_split_t bitv,
                                         unsigned int tip_count);

struct split_node_pair {
  pll_split_t split;
  pll_unode_t * node;
};

struct cb_split_params
{
  struct split_node_pair * split_nodes;
  unsigned int tip_count;
  unsigned int split_size;
  unsigned int split_len;
  unsigned int split_count;  /* number of splits already set */
  int *id_to_split;          /* map between node/subnode ids and splits */
};

/**
 * Check whether tip node indices in 2 trees are consistent to each other.
 *
 * Data pointers must contain the node index at first position
 * @param  t1 first tree
 * @param  t2 second tree
 *
 * @return PLL_SUCCESS if node indices are consistent,
 *         PLL_FAILURE otherwise (check pll_errmsg for details)
 */
PLL_EXPORT int pllmod_utree_consistency_check(pll_utree_t * t1,
                                              pll_utree_t * t2)
{
  unsigned int i;
  unsigned int node_id;
  int retval = PLL_SUCCESS;
  pll_unode_t ** tipnodes = t1->nodes;
  char ** tipnames;
  unsigned int tip_count = t1->tip_count;

  if (tip_count != t2->tip_count)
  {
    pllmod_set_error(PLLMOD_TREE_ERROR_INVALID_TREE,
                     "Trees do not have the same number of tips\n");
    return PLL_FAILURE;
  }


  tipnames = (char **) malloc (tip_count * sizeof(char *));
  if (!(tipnodes && tipnames))
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for tipnodes and tipnames\n");
    return PLL_FAILURE;
  }

  /* fill names table */
  for (i = 0; i < tip_count; ++i)
  {
    node_id = tipnodes[i]->node_index;
    tipnames[node_id] = tipnodes[i]->label;
  }

  /* check names consistency */
  tipnodes = t2->nodes;
  for (i = 0; i < tip_count; ++i)
  {
    node_id = tipnodes[i]->node_index;
    if (strcmp(tipnames[node_id], tipnodes[i]->label))
    {
      retval = PLL_FAILURE;
      break;
    }
  }

  free(tipnames);
  return retval;
}

/**
 * Set t2 tip node indices consistent with t1.
 *
 * Data pointers must contain the node index at first position
 * @param  t1 reference tree
 * @param  t2 second tree
 *
 * @return PLL_SUCCESS if the operation was applied correctly,
 *         PLL_FAILURE otherwise (check pll_errmsg for details)
 */
PLL_EXPORT int pllmod_utree_consistency_set(pll_utree_t * t1,
                                            pll_utree_t * t2)
{
  unsigned int i, j;
  unsigned int node_id;
  int retval = PLL_SUCCESS, checkval;
  pll_unode_t ** tipnodes = t1->nodes;
  char ** tipnames;
  unsigned int tip_count = t1->tip_count;

  if (tip_count != t2->tip_count)
  {
    pllmod_set_error(PLLMOD_TREE_ERROR_INVALID_TREE,
                     "Trees do not have the same number of tips\n");
    return PLL_FAILURE;
  }

  tipnames = (char **) malloc (tip_count * sizeof(char *));

  if (!(tipnodes && tipnames))
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for tipnodes and tipnames\n");
    return PLL_FAILURE;
  }

  /* fill names table */
  for (i = 0; i < tip_count; ++i)
  {
    node_id = tipnodes[i]->node_index;
    tipnames[node_id] = tipnodes[i]->label;
  }

  /* set names consistency */
  tipnodes = t2->nodes;
  for (i = 0; i < tip_count; ++i)
  {
    // node_id = get_utree_node_id(tipnodes[i]);
    pll_unode_t * tipnode = tipnodes[i];
    checkval = 0;
    for (j = 0; j < tip_count; ++j)
    {
      if (!strcmp(tipnames[j], tipnode->label))
      {
        checkval = 1;
        tipnode->node_index = j;
        break;
      }
    }
    if (!checkval)
    {
      retval = PLL_FAILURE;
      break;
    }
  }

  free(tipnames);
  return retval;
}

/******************************************************************************/
/* discrete operations */

PLL_EXPORT unsigned int pllmod_utree_rf_distance(pll_unode_t * t1,
                                                 pll_unode_t * t2,
                                                 unsigned int tip_count)
{
  unsigned int rf_distance;

  /* reset pll_error */
  pll_errno = 0;

  /* split both trees */
  pll_split_t * s1 = pllmod_utree_split_create(t1, tip_count, NULL);
  pll_split_t * s2 = pllmod_utree_split_create(t2, tip_count, NULL);

  /* compute distance */
  rf_distance = pllmod_utree_split_rf_distance(s1, s2, tip_count);

  /* clean up */
  pllmod_utree_split_destroy(s1);
  pllmod_utree_split_destroy(s2);

  assert(rf_distance <= 2*(tip_count-3));
  return rf_distance;
}

/*
 * Precondition: splits must be normalized and sorted!
 */
PLL_EXPORT unsigned int pllmod_utree_split_rf_distance(pll_split_t * s1,
                                                       pll_split_t * s2,
                                                       unsigned int tip_count)
{
  unsigned int split_count = tip_count - 3;
  unsigned int split_len   = bitv_length(tip_count);
  unsigned int equal = 0;
  unsigned int s1_idx = 0,
               s2_idx = 0;

  for (s1_idx=0; s1_idx<split_count && s2_idx<split_count; ++s1_idx)
  {
    int cmp = compare_splits(s1[s1_idx], s2[s2_idx], split_len);
    if (!cmp)
    {
      equal++;
      s2_idx++;
    }
    else
    {
      if (cmp > 0)
      {
        while(++s2_idx < split_count &&
              (cmp = compare_splits(s1[s1_idx], s2[s2_idx], split_len)) > 0);
        if (!cmp)
        {
           equal++;
           //s2_idx++;
        }
      }
    }
  }

  assert(equal <= (tip_count-3));

  return 2*(tip_count - 3 - equal);
}



/******************************************************************************/
/* tree split functions */

/*
 * Normalize and sort.
 * Warning! first position might change, so if splits were allocated together
 * you should keep a pointer to the original first position such that you can
 * deallocate it afterwards!
 */
/**
 * normalizes and sorts a set of splits
 *
 * @param s           set of splits
 * @param tip_count   number of tips
 * @param split_count numer of splits in 's'
 * @param keep_first  do not change first pointer in 's' (i.e., s[0])
 *
 * `keep_fist` parameter is important if the set of splits were allocated in
 * a contiguous chunk of memory and you want to use s[0] to deallocate it in
 * the future.
 */
PLL_EXPORT void pllmod_utree_split_normalize_and_sort(pll_split_t * s,
                                                      unsigned int tip_count,
                                                      unsigned int split_count,
                                                      int keep_first)
{
  unsigned int i;
  unsigned int split_len;

  pllmod_reset_error();

  pll_split_t first_split;
  for (i=0; i<split_count;++i)
    bitv_normalize(s[i], tip_count);

  first_split = s[0];
  qsort(s, split_count, sizeof(pll_split_t), _cmp_splits);

  if (keep_first && first_split != s[0])
  {
    split_len = bitv_length(tip_count);

    /* find first split */
    for (i=1; s[i] != first_split && i < split_count; ++i);
    assert(i < split_count);

    /* swap */
    void * aux_mem = malloc(sizeof(pll_split_base_t) * split_len);
    if (!aux_mem)
    {
      pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                       "Cannot allocate memory for auxiliary array\n");
      return;
    }

    memcpy(aux_mem, first_split,  sizeof(pll_split_base_t) * split_len);
    memcpy(first_split, s[0],     sizeof(pll_split_base_t) * split_len);
    memcpy(s[0], aux_mem, sizeof(pll_split_base_t) * split_len);
    free(aux_mem);
    s[i] = s[0];
    s[0] = first_split;
  }
}

PLL_EXPORT void pllmod_utree_split_show(pll_split_t split, unsigned int tip_count)
{
  unsigned int split_size   = sizeof(pll_split_base_t) * 8;
  unsigned int split_offset = tip_count % split_size;
  unsigned int split_len    = bitv_length(tip_count);
  unsigned int i, j;

  if (!split_offset) split_offset = split_size;

  for (i=0; i<(split_len-1); ++i)
    for (j=0; j<split_size; ++j)
      (split[i]&(1u<<j))?putchar('*'):putchar('-');
  for (j=0; j<split_offset; ++j)
    (split[i]&(1u<<j))?putchar('*'):putchar('-');
}

/*
 * Note: This function returns the splits according to the node indices at the tips!
 *
 * split_to_node_map can be NULL
 */
PLL_EXPORT pll_split_t * pllmod_utree_split_create(pll_unode_t * tree,
                                                   unsigned int tip_count,
                                                   pll_unode_t ** split_to_node_map)
{
  unsigned int i;
  unsigned int split_count, split_len, split_size;
  pll_split_t * split_list;   /* array with ordered split pointers */
  pll_split_t splits;         /* contiguous array of splits, as size is known */
  struct split_node_pair * split_nodes;
  pll_split_t first_split;

  /* as many non-trivial splits as inner branches */
  split_count   = tip_count - 3;
  split_size = (sizeof(pll_split_base_t) * 8);
  split_len  = (tip_count / split_size) + (tip_count % (sizeof(pll_split_base_t) * 8) > 0);

  split_list = (pll_split_t *) malloc(split_count * sizeof(pll_split_t));
  if (!split_list)
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for split list\n");
    return NULL;
  }

  split_nodes = (struct split_node_pair *) malloc(split_count * sizeof(struct split_node_pair));
  if (!split_nodes)
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for split-node pairs\n");
    free (split_list);
    return NULL;
  }

  splits = (pll_split_t) calloc(split_count * split_len, sizeof(pll_split_base_t));
  if (!splits)
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for splits\n");
    free (split_list);
    free (split_nodes);
    return NULL;
  }

  for (i=0; i<split_count; ++i)
  {
    split_nodes[i].split = splits + i*split_len;
    split_list[i] = splits + i*split_len;
  }

  struct cb_split_params split_data;
  split_data.split_nodes = split_nodes;
  split_data.split_len   = split_len;
  split_data.split_size  = split_size;
  split_data.tip_count   = tip_count;
  split_data.split_count = 0;

  /* reserve positions for node and subnode ids */
  split_data.id_to_split = (int *) malloc(sizeof(int) * 3 * (tip_count - 2));

  if (!split_data.id_to_split)
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for splits\n");
    return NULL;
  }

  for (i=0; i<3*(tip_count-2);++i)
    split_data.id_to_split[i] = -1;

  if (pllmod_utree_is_tip(tree))
    tree = tree->back;

  /* traverse for computing the scripts */
  pllmod_utree_traverse_apply(tree,
                              NULL,
                              NULL,
                              &cb_get_splits,
                              &split_data);

  assert(split_data.split_count == split_count);

  for (i=0; i<split_count; ++i)
  {

  }

  free(split_data.id_to_split);

  for (i=0; i<split_count; ++i)
    bitv_normalize(split_list[i], tip_count);

  /* sort map and split list together */
  first_split = split_nodes[0].split;
  qsort(split_nodes, split_count, sizeof(struct split_node_pair), _cmp_split_node_pair);

  /* if first item has changed, swap them such that the array can be deallocated */
  if (first_split != split_nodes[0].split)
  {
    /* find first split */
    for (i=1; split_nodes[i].split != first_split && i < split_count; ++i);
    assert(i < split_count);

    /* swap */
    void * aux_mem = malloc(sizeof(pll_split_base_t) * split_len);
    if (!aux_mem)
    {
      pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                       "Cannot allocate memory for auxiliary array\n");
      return NULL;
    }

    memcpy(aux_mem, first_split,  sizeof(pll_split_base_t) * split_len);
    memcpy(first_split, split_nodes[0].split,     sizeof(pll_split_base_t) * split_len);
    memcpy(split_nodes[0].split, aux_mem, sizeof(pll_split_base_t) * split_len);
    free(aux_mem);
    split_nodes[i].split = split_nodes[0].split;
    split_nodes[0].split = first_split;
  }

  for (i=0; i<split_count; ++i)
  {
    split_list[i] = split_nodes[i].split;
    assert(split_is_valid_and_normalized(split_list[i], tip_count));
  }

  /* update output arrays */
  if (split_to_node_map)
  {
    for (i=0; i<split_count; ++i)
    {
      split_list[i] = split_nodes[i].split;
      split_to_node_map[i] = split_nodes[i].node;
    }
  }
  else
  {
    for (i=0; i<split_count; ++i)
      split_list[i] = split_nodes[i].split;
  }

  free(split_nodes);

  return split_list;
}

/**
 * Creates or updates hashtable with splits (and their support)
 *
 * @param splits_hash    hashtable to update, NULL: create new hashtable
 * @param tip_count      number of tips
 * @param split_count    number of splits in 'splits'
 * @param support        support values for the split
 * @param update_only    0: insert new values as needed,
 *                       1: only increment support for existing splits
 *
 * @returns hashtable with splits
 */
PLL_EXPORT bitv_hashtable_t *
pllmod_utree_split_hashtable_insert(bitv_hashtable_t * splits_hash,
                                    pll_split_t * splits,
                                    unsigned int tip_count,
                                    unsigned int split_count,
                                    const double * support,
                                    int update_only)
{
  unsigned int i;

  if (!splits_hash)
  {
    /* create new hashtable */
    splits_hash = hash_init(tip_count * 10, tip_count);
    /* hashtable is empty, so update_only doesn't make sense here */
    update_only = 0;
  }

  if (!splits_hash)
  {
    return PLL_FAILURE;
  }

  /* insert splits */
  for (i=0; i<split_count; ++i)
  {
    if (update_only)
    {
      hash_update(splits[i],
                  splits_hash,
                  HASH_KEY_UNDEF,
                  support ? support[i] : 1.0,
                  0);
    }
    else
    {
      hash_insert(splits[i],
                  splits_hash,
                  i,
                  HASH_KEY_UNDEF,
                  support ? support[i] : 1.0,
                  0);
    }
  }

  return splits_hash;
}

PLL_EXPORT bitv_hash_entry_t *
pllmod_utree_split_hashtable_lookup(bitv_hashtable_t * splits_hash,
                                    pll_split_t split,
                                    unsigned int tip_count)
{
  unsigned int split_len = bitv_length(tip_count);
  hash_key_t position = hash_get_key(split, split_len) % splits_hash->table_size;
  bitv_hash_entry_t *p = splits_hash->table[position];

  for(; p!= NULL; p = p->next)
  {
    if(!compare_splits (p->bit_vector,
                        split,
                        split_len))
      return p;
  }

  return 0;
}

PLL_EXPORT
void pllmod_utree_split_hashtable_destroy(bitv_hashtable_t * hash)
{
  if (hash)
    hash_destroy(hash);
}

PLL_EXPORT void pllmod_utree_split_destroy(pll_split_t * split_list)
{
  free(split_list[0]);
  free(split_list);
}

/******************************************************************************/
/* static functions */

static inline void merge_split(pll_split_t to,
                         const pll_split_t from,
                         unsigned int split_len)
{
  unsigned int i;
  for (i=0;i<split_len;++i)
    to[i] |= from[i];
}



/**
 * Callback function for computing the splits at each branch
 * The splits will be stored in data->splits
 * at positions given by node index
 */
static int cb_get_splits(pll_unode_t * node, void *data)
{
  struct cb_split_params * split_data = (struct cb_split_params *) data;
  pll_split_t current_split;

  unsigned int tip_count     = split_data->tip_count;
  unsigned int split_size    = split_data->split_size;
  unsigned int split_len     = split_data->split_len;
  unsigned int my_split_id, child_split_id;
  unsigned int my_map_id, back_map_id;
  unsigned int tip_id, split_id;

  if (!(pllmod_utree_is_tip(node) || pllmod_utree_is_tip(node->back)))
  {
    my_map_id   = get_utree_splitmap_id(node, tip_count);
    back_map_id = get_utree_splitmap_id(node->back, tip_count);
    my_split_id = split_data->split_count;

    /* check if the split for the branch was already set */
    /* note that tree traversals visit the virtual root branch twice */
    if (split_data->id_to_split[my_map_id] >= 0)
    {
      return 1;
    }

    assert(my_split_id < (tip_count - 3));
    split_data->id_to_split[my_map_id] = (int) my_split_id;
    split_data->id_to_split[back_map_id] = (int) my_split_id;

    split_data->split_nodes[my_split_id].node = node;

    /* get current split to fill */
    current_split = split_data->split_nodes[my_split_id].split;
    /* increase number of splits */
    split_data->split_count++;

    /* add the split from left branch */
    if (!pllmod_utree_is_tip(node->next->back))
    {
      child_split_id = (unsigned int)
        split_data->id_to_split[get_utree_splitmap_id(node->next, tip_count)];

      memcpy(current_split, split_data->split_nodes[child_split_id].split,
             sizeof(pll_split_base_t) * split_len);
    }
    else
    {
      tip_id     = node->next->back->node_index;
      assert(tip_id < tip_count);
      split_id   = tip_id / split_size;
      tip_id    %= split_size;
      current_split[split_id] = (1 << tip_id);
    }

    /* add the split from right branch */
    if (!pllmod_utree_is_tip(node->next->next->back))
    {
      child_split_id = (unsigned int)
        split_data->id_to_split[get_utree_splitmap_id(
            node->next->next, tip_count)];
      merge_split(current_split, split_data->split_nodes[child_split_id].split, split_len);
    }
    else
    {
      tip_id     = node->next->next->back->node_index;
      assert(tip_id < tip_count);
      split_id   = tip_id / split_size;
      tip_id    %= split_size;
      current_split[split_id] |= (1 << tip_id);
    }
  }

  /* continue */
  return 1;
}

/*
 * The order of the splits is not really significant, as long as the two
 * following agree.
 *
 * _cmp_splits is used for sorting.
 * compare_splits is used for comparing splits from different trees
 */
static int compare_splits (pll_split_t s1,
                           pll_split_t s2,
                           unsigned int split_len)
{
  unsigned int i;

  for (i=0; i<split_len; ++i)
  {
    if (s1[i] != s2[i])
      return (int) (s1[i] > s2[i]?1:-1);
  }
  return 0;
}

/*
 * Precondition: splits *must* be different.
 */
static int _cmp_splits (const void * a, const void * b)
{
  const pll_split_t * s1 = (const pll_split_t *) a;
  const pll_split_t * s2 = (const pll_split_t *) b;
  unsigned int limit = 10000; /* max_taxa = split_size * 10^4 */
  int i = 0;
  for(;((*s1)[i]==(*s2)[i]) && limit;--limit,++i);
  assert(limit);
  return (int) ((*s1)[i]>(*s2)[i]?1:-1);
}

static int _cmp_split_node_pair (const void * a, const void * b)
{
  const struct split_node_pair * s1 = (const struct split_node_pair *) a;
  const struct split_node_pair * s2 = (const struct split_node_pair *) b;

  return (_cmp_splits(&s1->split, &s2->split));
}
/*
  The position of the node in the map of branches to splits is computed
  according to the node id.
 */
static unsigned int get_utree_splitmap_id(pll_unode_t * node,
                                          unsigned int tip_count)
{
  unsigned int node_id = node->node_index;
  assert(node_id >= tip_count);
  return node_id - tip_count;
}


/*
 * Returns 1 if the split is valid (not all 0s or all 1s) and normalized;
 * returns 0 otherwise
 * */
static int split_is_valid_and_normalized(const pll_split_t bitv,
                                         unsigned int tip_count)
{
  // this will also automatically check for all-0s case
  if (!bitv_is_normalized(bitv))
    return 0;

  // now check that we don't have all 1s
  unsigned int split_size  = sizeof(pll_split_base_t) * 8;
  unsigned int split_offset = tip_count % split_size;
  unsigned int split_len    = bitv_length(tip_count);
  unsigned int i = 0;
  unsigned int all1 = ~0;
  unsigned int mask = all1;
  for (i=0; i<split_len-1; ++i)
  {
    mask &= bitv[i];
  }
  if (split_offset)
    mask &= bitv[split_len-1] | ((1<<split_offset) - 1);
  else
    mask &= bitv[split_len-1];

  return (mask == all1) ? 0 : 1;
}
