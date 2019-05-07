/*
    Copyright (C) 2019 Sarah Lutteropp, Alexey Kozlov

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

    Contact: Sarah Lutteropp <Sarah.Lutteropp@h-its.org>,
    Heidelberg Institute for Theoretical Studies,
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "pll_tree.h"
#include "tree_hashtable.h"

#include "../pllmod_common.h"

typedef struct index_information
{
  unsigned int idx;
  unsigned int idx_left;
  unsigned int idx_right;
} index_information_t;

typedef struct tbe_data
{
  unsigned int* subtree_size;
  index_information_t* idx_infos;
  unsigned int* count_ones;
  unsigned int nodes_count;
  unsigned int trav_size;
  unsigned int tip_count;
  unsigned int tip_count_div_2;
} tbe_data_t;

int cb_full_traversal(pll_unode_t * node)
{
  (void) node;
  return 1;
}

void postorder_init_recursive(pll_unode_t * node, unsigned int * index,
                              unsigned int * subtree_size, index_information_t* idx_infos)
{
  if (node->next == NULL)
  {
    subtree_size[node->clv_index] = 1;
    return;
  }
  pll_unode_t * snode = node->next;
  do
  {
    postorder_init_recursive(snode->back, index, subtree_size, idx_infos);
    snode = snode->next;
  } while (snode && snode != node);
  index_information_t info;
  info.idx = node->clv_index;
  info.idx_left = node->next->back->clv_index;
  info.idx_right = node->next->next->back->clv_index;
  subtree_size[node->clv_index] = subtree_size[info.idx_left] + subtree_size[info.idx_right];
  idx_infos[*index] = info;
  *index = *index + 1;
}

void postorder_init(pll_unode_t * root, unsigned int * trav_size,
                    unsigned int * subtree_size, index_information_t* idx_infos)
{
  *trav_size = 0;
  postorder_init_recursive(root->back, trav_size, subtree_size, idx_infos);
  postorder_init_recursive(root, trav_size, subtree_size, idx_infos);
}

tbe_data_t* init_tbe_data(pll_unode_t * root, unsigned int tip_count)
{
  tbe_data_t* data = (tbe_data_t*) malloc(sizeof(tbe_data_t));
  data->tip_count = tip_count;
  data->tip_count_div_2 = tip_count / 2;
  data->trav_size = 0;
  data->nodes_count = 2 * tip_count - 2;
  data->subtree_size = (unsigned int*) malloc(sizeof(unsigned int) * data->nodes_count);
  data->idx_infos = (index_information_t*) malloc(sizeof(index_information_t) * data->nodes_count);
  data->count_ones = (unsigned int*) malloc(sizeof(unsigned int) * data->nodes_count);
  postorder_init(root, &data->trav_size, data->subtree_size, data->idx_infos);
  return data;
}

void free_tbe_data(tbe_data_t* data)
{
  free(data->subtree_size);
  free(data->idx_infos);
  free(data->count_ones);
  free(data);
}

unsigned int search_mindist(const pllmod_tbe_split_info_t* query,
                            tbe_data_t* data)
{
  unsigned int min_dist = query->p - 1;
  unsigned int* count_ones = data->count_ones;

  // initialize the leaf node informations...
  for (size_t i = 0; i < query->left_leaf_idx; ++i)
  {
    count_ones[i] = !query->subtree_res;
  }
  for (size_t i = query->left_leaf_idx; i <= query->right_leaf_idx; ++i)
  {
    count_ones[i] = query->subtree_res;
  }
  for (size_t i = query->right_leaf_idx + 1; i < data->nodes_count; ++i)
  {
    count_ones[i] = !query->subtree_res;
  }

  for (size_t i = 0; i < data->trav_size; ++i)
  {
    unsigned int idx = data->idx_infos[i].idx;
    unsigned int idx_left = data->idx_infos[i].idx_left;
    unsigned int idx_right = data->idx_infos[i].idx_right;
    count_ones[idx] = count_ones[idx_left] + count_ones[idx_right];
    unsigned int count_zeros = data->subtree_size[idx] - count_ones[idx];
    unsigned int dist_cand = query->p - count_zeros + count_ones[idx];

    if (dist_cand > data->tip_count_div_2)
    {
      dist_cand = data->tip_count - dist_cand;
    }
    if (dist_cand < min_dist)
    {
      min_dist = dist_cand;
      if (min_dist == 1)
      {
        return min_dist;
      }
    }
  }
  return min_dist;
}

/* This function computes a lower bound of Hamming distance between splits:
 * the computation terminates as soon as current distance values exceeds min_hdist.
 * This allows for substantial time savings if we are looking for the
 * minimum Hamming distance (as in TBE computation below).
 *
 * WARNING: This function does not check that hdist < N/2. Therefore,
 * it should be called twice, with original and inverted s1 (or s2),
 * to account for possible complementary split encoding.
 * */
static unsigned int utree_split_hamming_distance_lbound(pll_split_t s1,
                                                        pll_split_t s2,
                                                        unsigned int split_len,
                                                        unsigned int min_hdist)
{
  unsigned int hdist = 0;
  unsigned int i;

  for (i = 0; (i < split_len) && (hdist <= min_hdist); ++i)
  {
    hdist += PLL_POPCNT32(s1[i] ^ s2[i]);
  }

  return hdist;
}


/*
 *
 * API functions
 *
 */

PLL_EXPORT
pllmod_tbe_split_info_t * pllmod_utree_tbe_nature_init(pll_unode_t * ref_root,
                                                       unsigned int tip_count,
                                                       const pll_unode_t** split_to_node_map)
{
  unsigned int nodes_count = 2 * tip_count - 2;
  unsigned int split_count = tip_count - 3;
  unsigned int subtree_size[nodes_count];
  unsigned int a_leaf_idx[nodes_count];
  unsigned int b_leaf_idx[nodes_count];

  pllmod_tbe_split_info_t* split_info = NULL;
  pll_unode_t ** travbuffer = NULL;

  split_info = (pllmod_tbe_split_info_t*) malloc(sizeof(pllmod_tbe_split_info_t) * split_count);

  travbuffer = (pll_unode_t **) malloc(nodes_count * sizeof(pll_unode_t *));

  if (!split_info || !travbuffer)
  {
    if (split_info)
      free(split_info);
    if (travbuffer)
      free(travbuffer);
    pllmod_set_error(PLL_ERROR_MEM_ALLOC, "Cannot allocate memory\n");
    return PLL_FAILURE;
  }

  // do a post-order traversal of the reference tree.
  unsigned int trav_size;
  pll_utree_traverse(ref_root, PLL_TREE_TRAVERSE_POSTORDER, cb_full_traversal, travbuffer, &trav_size);
  for (unsigned int i = 0; i < trav_size; ++i)
  { // first, we compute the subtree sizes.
    unsigned int idx = travbuffer[i]->clv_index;
    if (travbuffer[i]->next == NULL)
    { // we are at a leaf node
      subtree_size[idx] = 1;
      a_leaf_idx[idx] = idx;
      b_leaf_idx[idx] = idx;
    }
    else
    {
      unsigned int idxLeft = travbuffer[i]->next->back->clv_index;
      unsigned int idxRight = travbuffer[i]->next->next->back->clv_index;
      subtree_size[idx] = subtree_size[idxLeft] + subtree_size[idxRight];
      a_leaf_idx[idx] = a_leaf_idx[idxLeft];
      b_leaf_idx[idx] = b_leaf_idx[idxRight];
    }
  }

  // now we need to check for each split in the reference tree if we need to do p=subtree_size or p = n-subtree_size
  for (unsigned int i = 0; i < split_count; ++i)
  {
    unsigned int node_idx = split_to_node_map[i]->clv_index;
    unsigned int first_leaf_idx = a_leaf_idx[node_idx];
    bool first_leaf_one = (first_leaf_idx > 0); // because taxon 0 is by convention always zero
    if (subtree_size[node_idx] <= tip_count - subtree_size[node_idx])
    {
      split_info[i].p = subtree_size[node_idx];
      split_info[i].subtree_res = !first_leaf_one ? first_leaf_one : !first_leaf_one;
    }
    else
    {
      split_info[i].p = tip_count - subtree_size[node_idx];
      split_info[i].subtree_res = first_leaf_one ? first_leaf_one : !first_leaf_one;
    }
    split_info[i].left_leaf_idx = a_leaf_idx[node_idx];
    split_info[i].right_leaf_idx = b_leaf_idx[node_idx];
  }
  free(travbuffer);

  return split_info;
}

PLL_EXPORT int pllmod_utree_tbe_nature(pll_split_t * ref_splits,
                                       pll_split_t * bs_splits,
                                       pll_unode_t* bs_root,
                                       unsigned int tip_count,
                                       double * support,
                                       pllmod_tbe_split_info_t* split_info)
{
  unsigned int i;
  unsigned int split_count = tip_count - 3;

  if (!ref_splits || !bs_splits || !support || !split_info)
  {
    pllmod_set_error(PLL_ERROR_PARAM_INVALID, "Parameter is NULL!\n");
    return PLL_FAILURE;
  }

  bitv_hashtable_t * bs_splits_hash = pllmod_utree_split_hashtable_insert(NULL, bs_splits,
                                                                          tip_count, split_count,
                                                                          NULL, 0);

  if (!bs_splits_hash)
    return PLL_FAILURE;

  tbe_data_t* tbe_data = NULL;

  /* iterate over all splits of the reference tree */
  for (i = 0; i < split_count; i++)
  {
    pll_split_t ref_split = ref_splits[i];

    if (pllmod_utree_split_hashtable_lookup(bs_splits_hash, ref_split, tip_count))
    {
      /* found identical split in a bootstrap tree -> assign full support */
      support[i] = 1.0;
      continue;
    }

    if (split_info[i].p == 2)
    { // no need for further searching
      support[i] = 0.0;
      continue;
    }

    if (!tbe_data)
      tbe_data = init_tbe_data(bs_root, tip_count);

    // else, we are in the search for minimum distance...
    unsigned int min_hdist = search_mindist(&split_info[i], tbe_data);
    //assert(min_hdist > 0);
    support[i] = 1.0 - (((double) min_hdist) / (split_info[i].p - 1));
  }

  pllmod_utree_split_hashtable_destroy(bs_splits_hash);

  if (tbe_data)
    free_tbe_data(tbe_data);

  return PLL_SUCCESS;
}


/* This is an old, naive and rather inefficient TBE computation method by Alexey,
 * keep it here just in case */
PLL_EXPORT int pllmod_utree_tbe_naive(pll_split_t * ref_splits,
                                      pll_split_t * bs_splits,
                                      unsigned int tip_count,
                                      double * support)
{
  unsigned int i, j, k;
  unsigned int split_count = tip_count - 3;
  unsigned int split_len = bitv_length(tip_count);
  unsigned int split_size  = sizeof(pll_split_base_t) * 8;
  unsigned int split_offset = tip_count % split_size;
  unsigned int split_mask = split_offset ? (1u<<split_offset) - 1 : ~0u;

  pll_split_t inv_split = NULL;
  unsigned int * bs_light = NULL;

  if (!ref_splits || !bs_splits || !support)
  {
    pllmod_set_error(PLL_ERROR_PARAM_INVALID, "Parameter is NULL!\n");
    return PLL_FAILURE;
  }

  bitv_hashtable_t * bs_splits_hash =
                          pllmod_utree_split_hashtable_insert(NULL,
                                                              bs_splits,
                                                              tip_count,
                                                              split_count,
                                                              NULL,
                                                              0);

  if (!bs_splits_hash)
  {
    return PLL_FAILURE;
  }

  inv_split = (pll_split_t) calloc(split_len, sizeof(pll_split_base_t));
  bs_light = calloc(split_count, sizeof(unsigned int));

  if (!inv_split || !bs_light)
  {
    if (inv_split)
      free(inv_split);
    if (bs_light)
      free(bs_light);
    pllmod_utree_split_hashtable_destroy(bs_splits_hash);
    pllmod_set_error(PLL_ERROR_MEM_ALLOC, "Cannot allocate memory\n");
    return PLL_FAILURE;
  }

  /* precompute lightside size for all bootstrap splits */
  for (j = 0; j < split_count; j++)
  {
    bs_light[j] = pllmod_utree_split_lightside(bs_splits[j], tip_count);
  }

  /* iterate over all splits of the reference tree */
  for (i = 0; i < split_count; i++)
  {
    pll_split_t ref_split = ref_splits[i];
    unsigned int p =  pllmod_utree_split_lightside(ref_split, tip_count);
    unsigned int min_hdist = p - 1;

    if (pllmod_utree_split_hashtable_lookup(bs_splits_hash, ref_split, tip_count))
    {
      /* found identical split in a bootstrap tree -> assign full support */
      support[i] = 1.0;
      continue;
    }

    /* inverse the reference split */
    for (k = 0; k < split_len; ++k)
    {
      inv_split[k] = ~ref_split[k];
    }
    /* clear unused bits in the last array element */
    inv_split[split_len-1] &= split_mask;

    /* iterate over all splits of the bootstrap tree */
    for (j = 0; j < split_count; j++)
    {
      unsigned int hdist, hdist_inv;

      /* this split is too far away -> skip it */
      if (abs(bs_light[j] - p) > min_hdist &&
          abs(tip_count - bs_light[j] - p) > min_hdist)
      {
        continue;
      }

      //      unsigned int hdist = pllmod_utree_split_hamming_distance(ref_split, bs_splits[j], tip_count);
      hdist = utree_split_hamming_distance_lbound(ref_split, bs_splits[j],
                                                  split_len, min_hdist);
      hdist_inv = utree_split_hamming_distance_lbound(inv_split, bs_splits[j],
                                                      split_len, min_hdist);
      min_hdist = PLL_MIN(min_hdist, hdist);
      min_hdist = PLL_MIN(min_hdist, hdist_inv);
    }

    assert(min_hdist > 0);

    support[i] = 1.0 - (((double) min_hdist) / (p - 1)) ;
  }

  pllmod_utree_split_hashtable_destroy(bs_splits_hash);
  free(inv_split);
  free(bs_light);

  return PLL_SUCCESS;
}

