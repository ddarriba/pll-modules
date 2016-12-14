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

#include "../pllmod_common.h"

static int cb_get_splits(pll_utree_t * node, void *data);
static inline void clone_split(pll_split_t to,
                         const pll_split_t from,
                         unsigned int split_len);
static inline void merge_split(pll_split_t to,
                         const pll_split_t from,
                         unsigned int split_len);
static void normalize_split(pll_split_t split, unsigned int tip_count);
static int _cmp_splits (const void * a, const void * b);
static int compare_splits (pll_split_t s1,
                           pll_split_t s2,
                           unsigned int split_len);
static unsigned int get_utree_splitmap_id(pll_utree_t * node,
                                          unsigned int tip_count);

struct cb_split_params
{
  pll_split_t * splits;
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
                                              pll_utree_t * t2,
                                              unsigned int tip_count)
{
  unsigned int i;
  unsigned int node_id;
  int retval = PLL_SUCCESS;
  pll_utree_t ** tipnodes;
  char ** tipnames;

  tipnodes = (pll_utree_t **) calloc ((size_t) tip_count, sizeof(pll_utree_t *));
  tipnames = (char **) malloc (tip_count * sizeof(char *));
  if (!(tipnodes && tipnames))
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for tipnodes and tipnames\n");
    return PLL_FAILURE;
  }

  pll_utree_query_tipnodes (t1, tipnodes);

  /* fill names table */
  for (i = 0; i < tip_count; ++i)
  {
    node_id = tipnodes[i]->node_index;
    tipnames[node_id] = tipnodes[i]->label;
  }

  /* check names consistency */
  pll_utree_query_tipnodes (t2, tipnodes);
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
  free(tipnodes);
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
                                            pll_utree_t * t2,
                                            unsigned int tip_count)
{
  unsigned int i, j;
  unsigned int node_id;
  int retval = PLL_SUCCESS, checkval;
  pll_utree_t ** tipnodes;
  char ** tipnames;

  tipnodes = (pll_utree_t **) calloc ((size_t) tip_count, sizeof(pll_utree_t *));
  tipnames = (char **) malloc (tip_count * sizeof(char *));

  if (!(tipnodes && tipnames))
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for tipnodes and tipnames\n");
    return PLL_FAILURE;
  }

  pll_utree_query_tipnodes (t1, tipnodes);

  /* fill names table */
  for (i = 0; i < tip_count; ++i)
  {
    node_id = tipnodes[i]->node_index;
    tipnames[node_id] = tipnodes[i]->label;
  }

  /* set names consistency */
  pll_utree_query_tipnodes (t2, tipnodes);
  for (i = 0; i < tip_count; ++i)
  {
    // node_id = get_utree_node_id(tipnodes[i]);
    pll_utree_t * tipnode = tipnodes[i];
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
  free(tipnodes);
  return retval;
}

/******************************************************************************/
/* discrete operations */

PLL_EXPORT unsigned int pllmod_utree_rf_distance(pll_utree_t * t1,
                                              pll_utree_t * t2,
                                              unsigned int tip_count)
{
  unsigned int split_count;
  unsigned int rf_distance;

  /* reset pll_error */
  pll_errno = 0;

  /* split both trees */
  pll_split_t * s1 = pllmod_utree_split_create(t1, tip_count, &split_count);
  pll_split_t * s2 = pllmod_utree_split_create(t2, tip_count, &split_count);

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
  unsigned int split_size  = (sizeof(pll_split_base_t) * 8);
  unsigned int split_len   = (tip_count / split_size) +
                            (tip_count % (sizeof(pll_split_base_t) * 8) > 0);
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
PLL_EXPORT void pllmod_utree_split_normalize_and_sort(pll_split_t * s,
                                                      unsigned int tip_count,
                                                      unsigned int split_count)
{
  unsigned int i;
  for (i=0; i<split_count;++i)
    normalize_split(s[i], tip_count);

  qsort(s, split_count, sizeof(pll_split_t), _cmp_splits);
}

PLL_EXPORT void pllmod_utree_split_show(pll_split_t split, unsigned int tip_count)
{
  unsigned int split_size = sizeof(pll_split_base_t) * 8;
  unsigned int split_offset = tip_count % split_size;
  unsigned int split_len  = tip_count / split_size + (split_offset>0);
  unsigned int i, j;

  for (i=0; i<(split_len-1); ++i)
    for (j=0; j<split_size; ++j)
      (split[i]&(1u<<j))?putchar('*'):putchar('-');
  for (j=0; j<split_offset; ++j)
    (split[i]&(1u<<j))?putchar('*'):putchar('-');
  printf("\n");
}

/*
 * Note: This function returns the splits according to the node indices at the tips!
 */
PLL_EXPORT pll_split_t * pllmod_utree_split_create(pll_utree_t * tree,
                                                   unsigned int tip_count,
                                                   unsigned int * _split_count)
{
  unsigned int i;
  unsigned int split_count, split_len, split_size;
  pll_split_t * split_list;
  pll_split_t splits;

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

  splits = (pll_split_t) calloc(split_count * split_len, sizeof(pll_split_base_t));

  if (!splits)
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for splits\n");
    free (split_list);
    return NULL;
  }

  for (i=0; i<split_count; ++i)
    split_list[i] = splits + i*split_len;

  struct cb_split_params split_data;
  split_data.split_len   = split_len;
  split_data.split_size  = split_size;
  split_data.splits      = split_list;
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
  free(split_data.id_to_split);

  /* normalize the splits such that first position is set */
  for (i=0; i<split_count;++i)
    normalize_split(split_list[i], tip_count);

  /* sort splits and keep first position pointing to the allocated array */
  pll_split_t aux = split_list[0];

  qsort(split_list, split_count, sizeof(pll_split_t), _cmp_splits);

#ifdef ULTRADEBUG
  /* check */
  for (i=1; i<split_count; ++i)
    if (compare_splits(split_list[i], split_list[i-1], split_len) <= 0)
    {
      printf("Error %d %d\n", i, _cmp_splits(&split_list[i], &split_list[i-1]));
      pllmod_utree_split_show(split_list[i], tip_count);
      pllmod_utree_split_show(split_list[i-1], tip_count);
    }
#endif

  for (i=0; split_list[i] != aux && i < split_count; ++i);
  assert(i < split_count);
  if (i>0)
  {
    /* swap */
    void * aux_mem = malloc(sizeof(pll_split_base_t) * split_len);
    if (!aux_mem)
    {
      pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                       "Cannot allocate memory for auxiliary array\n");
      return NULL;
    }

    memcpy(aux_mem, aux,           sizeof(pll_split_base_t) * split_len);
    memcpy(aux, split_list[0],     sizeof(pll_split_base_t) * split_len);
    memcpy(split_list[0], aux_mem, sizeof(pll_split_base_t) * split_len);
    free(aux_mem);
    split_list[i] = split_list[0];
    split_list[0] = aux;
  }

  if (_split_count) *_split_count = split_count;
  return split_list;
}

PLL_EXPORT void pllmod_utree_split_destroy(pll_split_t * split_list)
{
  free(split_list[0]);
  free(split_list);
}

/******************************************************************************/
/* static functions */

static inline void clone_split(pll_split_t to,
                         const pll_split_t from,
                         unsigned int split_len)
{
  memcpy(to, from, sizeof(pll_split_base_t) * split_len);
}

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
static int cb_get_splits(pll_utree_t * node, void *data)
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

    /* get current split to fill */
    current_split = split_data->splits[my_split_id];
    /* increase number of splits */
    split_data->split_count++;

    /* add the split from left branch */
    if (!pllmod_utree_is_tip(node->next->back))
    {
      child_split_id = (unsigned int)
        split_data->id_to_split[get_utree_splitmap_id(node->next, tip_count)];
      clone_split(current_split, split_data->splits[child_split_id], split_len);
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
      merge_split(current_split, split_data->splits[child_split_id], split_len);
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

static void normalize_split(pll_split_t split, unsigned int tip_count)
{
  unsigned int split_size = sizeof(pll_split_base_t) * 8;
  unsigned int split_offset = tip_count % split_size;
  unsigned int split_len  = tip_count / split_size + (split_offset>0);
  unsigned int i;

  int normalized = split[0]&1;

  if (!normalized)
  {
    for (i=0; i<split_len; ++i)
      split[i] = ~split[i];

    unsigned int mask = (1<<split_offset) - 1;
    split[split_len - 1] &= mask;
  }
}

/*
  The position of the node in the map of branches to splits is computed
  according to the node id.
 */
static unsigned int get_utree_splitmap_id(pll_utree_t * node,
                                          unsigned int tip_count)
{
  unsigned int node_id = node->node_index;
  assert(node_id >= tip_count);
  return node_id - tip_count;
}
