/*
 * tbe_functions.c
 *
 *  Created on: Feb 22, 2019
 *      Author: Sarah Lutteropp
 */

#include "pll_tree.h"

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
  unsigned int nodes_count;
  unsigned int trav_size;
  unsigned int tip_count;
  unsigned int tip_count_div_2;
} tbe_data_t;

int cb_full_traversal(pll_unode_t * node) {
  (void) node;
  return 1;
}

PLL_EXPORT refsplit_info_t* init_ref_splits_tbe(pll_unode_t * ref_root, unsigned int tip_count, const pll_unode_t** split_to_node_map) {
  unsigned int nodes_count = 2 * tip_count - 2;
  unsigned int split_count = tip_count - 3;
  refsplit_info_t* split_info = (refsplit_info_t*) malloc(sizeof(refsplit_info_t) * split_count);
  unsigned int subtree_size[nodes_count];
  unsigned int a_leaf_idx[nodes_count];
  unsigned int b_leaf_idx[nodes_count];
  // do a post-order traversal of the reference tree.
  pll_unode_t ** travbuffer = (pll_unode_t **) malloc(nodes_count * sizeof(pll_unode_t *));
  unsigned int trav_size;
  pll_utree_traverse(ref_root, PLL_TREE_TRAVERSE_POSTORDER, cb_full_traversal, travbuffer, &trav_size);
  for (unsigned int i = 0; i < trav_size; ++i)
  { // first, we compute the subtree sizes.
    unsigned int idx = travbuffer[i]->clv_index;
    if (travbuffer[i]->next == NULL) { // we are at a leaf node
      subtree_size[idx] = 1;
      a_leaf_idx[idx] = idx;
      b_leaf_idx[idx] = idx;
    } else {
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
    split_info[i].right_leaf_index = b_leaf_idx[node_idx];
  }
  free(travbuffer);

  return split_info;
}


void postorder_init_recursive(pll_unode_t * node, unsigned int * index, unsigned int * subtree_size, index_information_t* idx_infos) {
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

void postorder_init(pll_unode_t * root, unsigned int * trav_size, unsigned int * subtree_size, index_information_t* idx_infos) {
  *trav_size = 0;
  postorder_init_recursive(root->back, trav_size, subtree_size, idx_infos);
  postorder_init_recursive(root, trav_size, subtree_size, idx_infos);
}

tbe_data_t* init_tbe_data(pll_unode_t * root, unsigned int tip_count) {
  tbe_data_t* data = (tbe_data_t*) malloc(sizeof(tbe_data_t));
  data->tip_count = tip_count;
  data->tip_count_div_2 = tip_count / 2;
  data->trav_size = 0;
  data->nodes_count = 2 * tip_count - 2;
  data->subtree_size = (unsigned int*) malloc(sizeof(unsigned int) * data->nodes_count);
  data->idx_infos = (index_information_t*) malloc(sizeof(index_information_t) * data->nodes_count);
  postorder_init(root, &data->trav_size, data->subtree_size, data->idx_infos);
  return data;
}

void free_tbe_data(tbe_data_t* data) {
  free(data->subtree_size);
  free(data->idx_infos);
  free(data);
}

unsigned int search_mindist(const refsplit_info_t* query, const tbe_data_t* data, unsigned int* count_ones) {
  unsigned int min_dist = query->p - 1;
  // initialize the leaf node informations...
  for (size_t i = 0; i < query->left_leaf_idx; ++i)
  {
    count_ones[i] = !query->subtree_res;
  }
  for (size_t i = query->left_leaf_idx; i <= query->right_leaf_index; ++i)
  {
    count_ones[i] = query->subtree_res;
  }
  for (size_t i = query->right_leaf_index + 1; i < data->nodes_count; ++i)
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

unsigned int search_mindist_non_reusable_counts_array(const refsplit_info_t* query, const tbe_data_t* data) {
  unsigned int* count_ones = (unsigned int*) malloc(sizeof(unsigned int) * data->nodes_count);
  unsigned int res = search_mindist(query, data, count_ones);
  free(count_ones);
  return res;
}

PLL_EXPORT int pllmod_utree_split_transfer_support_nature(pll_split_t * ref_splits, pll_split_t * bs_splits, pll_unode_t* bs_root,
		unsigned int tip_count, double * support, refsplit_info_t* split_info) {
  unsigned int i;
  unsigned int split_count = tip_count - 3;

  if (!ref_splits || !bs_splits || !support || !split_info)
  {
    //pllmod_set_error(PLL_ERROR_PARAM_INVALID, "Parameter is NULL!\n");
    return PLL_FAILURE;
  }

  bitv_hashtable_t * bs_splits_hash = pllmod_utree_split_hashtable_insert(NULL, bs_splits, tip_count, split_count,
	NULL, 0);

  if (!bs_splits_hash)
  {
    return PLL_FAILURE;
  }

  tbe_data_t* tbe_data = NULL;
  unsigned int* count_ones = NULL;

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
    if (tbe_data == NULL)
    {
      tbe_data = init_tbe_data(bs_root, tip_count);
      count_ones = (unsigned int*) malloc(sizeof(unsigned int) * tbe_data->nodes_count);
    }
    // else, we are in the search for minimum distance...
    unsigned int min_hdist = search_mindist(&split_info[i], tbe_data, count_ones);
    //assert(min_hdist > 0);
    support[i] = 1.0 - (((double) min_hdist) / (split_info[i].p - 1));
  }

  pllmod_utree_split_hashtable_destroy(bs_splits_hash);

  if (tbe_data != NULL)
  {
    free_tbe_data(tbe_data);
    free(count_ones);
  }
  return PLL_SUCCESS;
}
