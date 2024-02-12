/*
 Copyright (C) 2021 Alexey Kozlov

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

 Contact: Alexey Kozlov <Alexey.Kozlov@h-its.org>,
 Exelixis Lab, Heidelberg Instutute for Theoretical Studies
 Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
 */

#include "pll_tree.h"
#include "tree_hashtable.h"

#include "../pllmod_common.h"

static inline unsigned int split_popcount(const pll_split_t bitv,
                                          unsigned int bit_count,
                                          unsigned int split_len)
{
  unsigned int split_size  = sizeof(pll_split_base_t) * 8;
  unsigned int setb = 0;
  unsigned int i;

  if (!split_len)
    split_len = bitv_length(bit_count);

  for (i = 0; i < split_len; ++i)
  {
    setb += (unsigned int) PLL_POPCNT32(bitv[i]);
  }

  /* IMPORTANT: correct for padding bits in the last element! */
  unsigned int split_offset = bit_count % split_size;
  if (split_offset)
  {
    unsigned int mask = (1<<split_offset) - 1;
    unsigned int last = bitv[split_len - 1];
    /* count set bits in the padding part of the bit vector */
    last &= ~mask;
    setb -= (unsigned int) PLL_POPCNT32(last);
  }

  return setb;
}

static inline void invert_split(pll_split_t bitv, unsigned int bit_count)
{
  unsigned int split_size  = sizeof(pll_split_base_t) * 8;
  unsigned int split_offset = bit_count % split_size;
  unsigned int split_len    = bitv_length(bit_count);
  unsigned int i;

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

static inline void copy_split(pll_split_t to,
                              const pll_split_t from,
                              unsigned int split_len)
{
  memcpy(to, from, split_len * sizeof(pll_split_base_t));
}

static inline void merge_split(pll_split_t to,
                         const pll_split_t from,
                         unsigned int split_len)
{
  unsigned int i;
  for (i=0;i<split_len;++i)
    to[i] |= from[i];
}

static inline int disjoint_split(const pll_split_t split1,
                                 const pll_split_t split2,
                                 unsigned int split_len)
{
  unsigned int i;
  for (i=0; i<split_len; ++i)
  {
    if (split1[i] & split2[i])
      return 0;
  }
  return 1;
}

static inline int empty_split(pll_split_t split,
                                   unsigned int split_len,
                                   unsigned int tip_count)
{
  unsigned int i;
  for (i=0;i<split_len-1;++i)
  {
    if (split[i])
      return 0;
  }

  unsigned int split_size  = sizeof(pll_split_base_t) * 8;
  unsigned int split_offset = tip_count % split_size;
  pll_split_base_t last_elem = split[split_len-1];
  if (split_offset)
  {
    unsigned int mask = (1<<split_offset) - 1;

    return ((last_elem & mask) == 0);
  }
  else
    return last_elem == 0;
}

static inline int full_split(pll_split_t split,
                                  unsigned int split_len,
                                  unsigned int tip_count)
{
  pll_split_base_t f = ~0;
  unsigned int i;
  for (i=0;i<split_len;++i)
  {
    if (split[i] != f)
      return 0;
  }

  unsigned int split_size  = sizeof(pll_split_base_t) * 8;
  unsigned int split_offset = tip_count % split_size;
  pll_split_base_t last_elem = split[split_len-1];
  if (split_offset)
  {
    unsigned int mask = (1<<split_offset) - 1;

    return ((last_elem & mask) == f);
  }
  else
    return last_elem == f;
}


static void truncate_splits(pll_split_set_t * split_set, unsigned int new_tip_count)
{
  unsigned int bitv_elem = split_set->split_size;
  unsigned int bitv_size = new_tip_count / bitv_elem;
  unsigned int bitv_off = new_tip_count % bitv_elem;
  if (bitv_off > 0)
    bitv_size++;

  if (split_set->tip_count > new_tip_count)
  {
    unsigned int mask = 0;
    for (unsigned int i = 0; i < (bitv_off ? bitv_off : bitv_elem); ++i)
      mask |= (1u << i);

    for (unsigned int i = 0; i < split_set->split_count; ++i)
    {
      pll_split_t last_elem = split_set->splits[i] + bitv_size - 1;
      *last_elem &= mask;
    }
  }
}

inline const pll_split_t get_node_split(const pll_split_t * splits,
                                        const pll_unode_t * node)
{
  return splits[node->node_index];
}

static const pll_split_t find_nonempty_regraft_split(pll_split_t * splits,
                                                     unsigned int split_len,
                                                     unsigned int tip_count,
                                                     const pll_split_t prune_split,
                                                     pll_unode_t * r_edge)
{
  pll_split_t regraft_split = NULL;
  pll_unode_t * left_node = r_edge;
  pll_unode_t * right_node = r_edge->back;

  while (1)
  {
    pll_split_t left_split = get_node_split(splits, left_node);
    pll_split_t right_split = get_node_split(splits, right_node);

    if (empty_split(right_split, split_len, tip_count) && !pllmod_utree_is_tip(left_node))
    {
      /* right subtree empty -> traverse left subtree */
      right_node = left_node->next->next->back;
      left_node = left_node->next->back;
    }
    else if (empty_split(left_split, split_len, tip_count) && !pllmod_utree_is_tip(right_node))
    {
      /* left subtree empty -> traverse right subtree */
      left_node = right_node->next->back;
      right_node = right_node->next->next->back;
    }
    else
    {
      /* non-trivial split found */
      if (disjoint_split(prune_split, right_split, split_len))
        regraft_split = right_split;
      else
        regraft_split = left_split;
      break;
    }
  }

  return regraft_split;
}

static int cb_get_all_splits(pll_unode_t * node, void *data)
{
  pll_split_t current_split, back_split;
  unsigned int my_split_id, back_split_id, child_split_id;
  unsigned int tip_id, split_id;

  pll_split_set_t * split_data = (pll_split_set_t *) data;
  unsigned int tip_count       = split_data->tip_count;
  unsigned int split_size      = split_data->split_size;
  unsigned int split_len       = split_data->split_len;

  my_split_id   = node->node_index;
  back_split_id = node->back->node_index;

  /* check if the split for the branch was already set */
  /* note that tree traversals visit the virtual root branch twice */
  if (split_data->id_to_split[my_split_id] >= 0)
    return 1;

  /* get current split to fill */
  current_split = split_data->splits[my_split_id];

  memset(current_split, 0, sizeof(pll_split_base_t) * split_len);

  if (pllmod_utree_is_tip(node))
  {
    /* trivial split */
    tip_id     = node->node_index;
    assert(tip_id < tip_count);
    split_id   = tip_id / split_size;
    tip_id    %= split_size;
    current_split[split_id] = (1 << tip_id);
  }
  else
  {
    /* add the split from branches */
    pll_unode_t * snode = node->next;
    while(snode != node)
    {
      child_split_id = snode->back->node_index;

      merge_split(current_split, split_data->splits[child_split_id], split_len);

      snode = snode->next;
    }
  }

  /* compute split in opposite direction -> simply invert all bits */
  back_split = split_data->splits[back_split_id];
  copy_split(back_split, current_split, split_len);
  invert_split(back_split, tip_count);

  /* here the mapping is trivial, we just use it to flag processed branches (see above) */
  split_data->id_to_split[my_split_id] = my_split_id;
  split_data->id_to_split[back_split_id] = back_split_id;

  /* increase number of splits -> two splits per branch! */
  split_data->split_count += 2;

  /* continue */
  return 1;
}

PLL_EXPORT pll_split_set_t * pllmod_utree_splitset_create(const pll_utree_t * tree)
{
  pll_split_set_t * split_set = (pll_split_set_t *) calloc(1, sizeof(pll_split_set_t));

  if (!split_set)
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for split list\n");
    return NULL;
  }

  /* init constraint */
  split_set->tip_count = tree->tip_count;
  split_set->split_size = sizeof(pll_split_base_t) * 8;
  split_set->split_len = bitv_length(tree->tip_count);
  split_set->split_count = tree->edge_count - tree->tip_count;
  split_set->splits = pllmod_utree_split_create(tree->vroot, tree->tip_count, NULL);

  if (!split_set->splits)
  {
    free(split_set);
    return NULL;
  }

  return split_set;
}

PLL_EXPORT pll_split_set_t * pllmod_utree_splitset_create_all(const pll_utree_t * tree)
{
  unsigned int i;
  pll_split_t split_storage;         /* contiguous array of splits, as size is known */

  pll_split_set_t * split_set = (pll_split_set_t *) calloc(1, sizeof(pll_split_set_t));

  if (!split_set)
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for split set\n");
    return NULL;
  }

  /* init constraint */
  split_set->tip_count = tree->tip_count;
  split_set->split_size = sizeof(pll_split_base_t) * 8;
  split_set->split_len = bitv_length(tree->tip_count);
  /* directed splits => two splits per branch */
  split_set->split_count = tree->edge_count * 2;

  split_set->splits = (pll_split_t *) malloc(split_set->split_count * sizeof(pll_split_t));
  if (!split_set->splits)
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for split list\n");
    pllmod_utree_splitset_destroy(split_set);
    return NULL;
  }

  split_storage = (pll_split_t) calloc(split_set->split_count * split_set->split_len,
                                       sizeof(pll_split_base_t));
  if (!split_storage)
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for splits\n");
    pllmod_utree_splitset_destroy(split_set);
    return NULL;
  }

  for (i = 0; i < split_set->split_count; ++i)
  {
    split_set->splits[i] = split_storage + i*split_set->split_len;
  }

  /* reserve positions for subnode ids */
  split_set->id_to_split = (int *) malloc(sizeof(int) * split_set->split_count);

  if (!split_set->id_to_split)
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for splits\n");
    pllmod_utree_splitset_destroy(split_set);
    return NULL;
  }

  pllmod_utree_splitset_update_all(split_set, tree);

  return split_set;
}

PLL_EXPORT int pllmod_utree_splitset_update_all(pll_split_set_t * split_set, const pll_utree_t * tree)
{
  unsigned int i;
  unsigned int tip_count = tree->tip_count;
  unsigned int split_count = tree->edge_count * 2;
  const pll_unode_t * vroot = pllmod_utree_is_tip(tree->vroot) ? tree->vroot->back : tree->vroot;

  if (split_set->tip_count != tip_count || split_set->split_count != split_count)
  {
    pllmod_set_error(PLL_ERROR_TREE_INVALID,
                     "Unable to update splits: tree has different number of tips/edges\n");
    return PLL_FAILURE;
  }

  /* clear branch processing flags */
  for (i = 0; i < split_set->split_count; ++i)
    split_set->id_to_split[i] = -1;

  split_set->split_count = 0;

  vroot = tree->vroot;
  if (pllmod_utree_is_tip(vroot))
    vroot = vroot->back;

  /* traverse for computing the scripts */
  pllmod_utree_traverse_apply((pll_unode_t *) vroot,
                              NULL,
                              NULL,
                              &cb_get_all_splits,
                              split_set);

  assert(split_set->split_count == split_count);

  return PLL_SUCCESS;
}

PLL_EXPORT void pllmod_utree_splitset_destroy(pll_split_set_t * split_set)
{
  if (split_set)
  {
    pllmod_utree_split_destroy(split_set->splits);
    free(split_set->id_to_split);
    free(split_set);
  }
}

PLL_EXPORT int pllmod_utree_constraint_check_spr(pll_split_set_t * cons_splits,
                                                 pll_split_set_t * tree_splits,
                                                 pll_unode_t * p_edge,
                                                 pll_unode_t * r_edge)
{
  int retval = PLL_SUCCESS;
  if (cons_splits)
  {
    unsigned int cons_tip_count = cons_splits->tip_count;
    unsigned int cons_split_len = cons_splits->split_len;
    pll_split_t * splits = tree_splits->splits;

    /* IDEA: check if the new branch added by SPR move contradicts the the topological constraint.
     * It gets a bit tricky if the branch immediately adjacent to the regrafting point is trivial
     * w.r.t. constraint, i.e. one of the subtrees does only contain a single constrained taxon.
     * In this case. we traverse into the larger subtree (which contains n-1 constraint taxa)
     * until we find the first non-trivial branch.
     */

    pll_split_t regraft_split = NULL;
    const pll_split_t prune_split = get_node_split(splits, p_edge->back);
    pll_split_t new_split = (pll_split_t) calloc(1, cons_split_len * sizeof(pll_split_base_t));

    unsigned int pruned_count = split_popcount(prune_split, cons_tip_count, cons_split_len);
    assert(pruned_count <= cons_tip_count);

    if (pruned_count < cons_tip_count-1)
    {
      /* remaining subtree contains at least 2 constrained taxa -> traverse into regraft subtree */
      regraft_split = find_nonempty_regraft_split(splits, cons_split_len, cons_tip_count, prune_split, r_edge);

      assert(!empty_split(regraft_split, cons_split_len, cons_tip_count));

      copy_split(new_split, regraft_split, cons_split_len);
      merge_split(new_split, prune_split, cons_split_len);
    }
    else
    {
      /* remaining subtree contains just 1 constrained taxon  -> traverse into pruned subtree */
      copy_split(new_split, prune_split, cons_split_len);
      invert_split(new_split, cons_tip_count);

      regraft_split = find_nonempty_regraft_split(splits, cons_split_len, cons_tip_count, new_split, p_edge);

      merge_split(new_split, regraft_split, cons_split_len);
    }

    /* check that newly introduced split is compatible with *all* constraint splits */
    for (unsigned int z = 0; z < cons_splits->split_count; ++z)
    {
       if (!pllmod_utree_compatible_splits(new_split, cons_splits->splits[z], cons_split_len, cons_tip_count))
       {
         retval = PLL_FAILURE;
         break;
       }
    }

    free(new_split);
  }

  return retval;
}


PLL_EXPORT int pllmod_utree_constraint_check_splits(pll_split_set_t * cons_splits, pll_split_set_t * tree_splits)
{
  int retval = PLL_SUCCESS;

  bitv_hashtable_t* splits_hash = pllmod_utree_split_hashtable_insert(NULL,
                                                                      tree_splits->splits,
                                                                      cons_splits->tip_count,
                                                                      tree_splits->split_count,
                                                                      NULL,
                                                                      0);

  for (size_t i = 0; i < cons_splits->split_count; ++i)
  {
    if (!pllmod_utree_split_hashtable_lookup(splits_hash, cons_splits->splits[i], cons_splits->tip_count))
    {
//      pllmod_utree_split_show(cons_splits->splits[i], cons_splits->tip_count);
//      printf("\n");
      retval = PLL_FAILURE;
      break;
    }
  }

  pllmod_utree_split_hashtable_destroy(splits_hash);

  return retval;
}

PLL_EXPORT int pllmod_utree_constraint_check_splits_tree(pll_split_set_t * cons_splits,
                                                         const pll_utree_t * tree)
{
  int retval = PLL_SUCCESS;

  pll_split_set_t * tree_splits =  pllmod_utree_splitset_create(tree);

  truncate_splits(tree_splits, cons_splits->tip_count);

  retval = pllmod_utree_constraint_check_splits(cons_splits, tree_splits);

  pllmod_utree_splitset_destroy(tree_splits);

  return retval;
}

PLL_EXPORT int pllmod_utree_constraint_check_tree(const pll_utree_t * cons_tree,
                                                  const pll_utree_t * tree)
{
  int retval = PLL_SUCCESS;

  pll_split_set_t * cons_splits =  pllmod_utree_splitset_create(cons_tree);

  retval = pllmod_utree_constraint_check_splits_tree(cons_splits, tree);

  pllmod_utree_splitset_destroy(cons_splits);

  return retval;
}

PLL_EXPORT int pllmod_utree_constraint_subtree_affected(const pll_split_set_t * cons_splits,
                                                        const pll_split_set_t * tree_splits,
                                                        pll_unode_t * p_edge)
{
  int retval = PLL_SUCCESS;

  const pll_split_t prune_split = get_node_split(tree_splits->splits, p_edge->back);
  unsigned int cons_split_len = cons_splits->split_len;
  unsigned int cons_tip_count = cons_splits->tip_count;

  /* zero constrained taxa in pruned OR remaining subtree */
  retval = (empty_split(prune_split, cons_split_len, cons_tip_count) ||
            full_split(prune_split, cons_split_len, cons_tip_count)) ? PLL_FAILURE : PLL_SUCCESS;

  return retval;
}
