#include "pll_tree.h"

static int cb_get_splits(pll_utree_t * node, void *data);
static inline void clone_split(pll_split_t to,
                         const pll_split_t from,
                         unsigned int split_len);
static inline void merge_split(pll_split_t to,
                         const pll_split_t from,
                         unsigned int split_len);
static void normalize_split(pll_split_t split, unsigned int n_tips);
static int _cmp_splits (const void * a, const void * b);
static int compare_splits (pll_split_t s1,
                           pll_split_t s2,
                           unsigned int split_len);

struct cb_split_params
{
  pll_split_t * splits;
  unsigned int tip_count;
  unsigned int split_size;
  unsigned int split_len;
};


/******************************************************************************/
/* discrete operations */

PLL_EXPORT unsigned int pll_utree_rf_distance(pll_utree_t * t1,
                                              pll_utree_t * t2,
                                              unsigned int n_tips)
{
  unsigned int n_splits;
  unsigned int rf_distance;

  /* split both trees */
  pll_split_t * s1 = pll_utree_split_create(t1, n_tips, &n_splits);
  pll_split_t * s2 = pll_utree_split_create(t2, n_tips, &n_splits);

  /* compute distance */
  rf_distance = pll_utree_rf_split_distance(s1, s2, n_tips);

  /* clean up */
  pll_utree_split_destroy(s1);
  pll_utree_split_destroy(s2);

  assert(rf_distance < 2*(n_tips-3));
  return rf_distance;
}

/*
 * Precondition: splits must be normalized and sorted!
 */
PLL_EXPORT unsigned int pll_utree_rf_split_distance(pll_split_t * s1,
                                                    pll_split_t * s2,
                                                    unsigned int n_tips)
{
  unsigned int n_splits = n_tips - 3;
  unsigned int split_size = (sizeof(pll_split_base_t) * 8);
  unsigned int split_len  = (n_tips / split_size) +
                            (n_tips % (sizeof(pll_split_base_t) * 8) > 0);
  unsigned int equal = 0;
  unsigned int s1_idx = 0,
               s2_idx = 0;

  for (s1_idx=0; s1_idx<n_splits && s2_idx<n_splits; ++s1_idx)
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
        while(++s2_idx < n_splits &&
              (cmp =compare_splits(s1[s1_idx], s2[s2_idx], split_len)) > 0);
        if (!cmp)
        {
           equal++;
           //s2_idx++;
        }
      }
    }
  }

  assert(equal <= (n_tips-3));

  return 2*(n_tips - 3 - equal);
}



/******************************************************************************/
/* tree split functions */

PLL_EXPORT void pll_utree_show_split(pll_split_t split, unsigned int n_tips)
{
  unsigned int split_size = sizeof(pll_split_base_t) * 8;
  unsigned int split_offset = n_tips % split_size;
  unsigned int split_len  = n_tips / split_size + (split_offset>0);
  unsigned int i, j;

  for (i=0; i<(split_len-1); ++i)
    for (j=0; j<split_size; ++j)
      (split[i]&(1u<<j))?putchar('*'):putchar('-');
  for (j=0; j<split_offset; ++j)
    (split[i]&(1u<<j))?putchar('*'):putchar('-');
  printf("\n");
}

/*
 * Note: This function returns the splits according to the p-matrix indices at the tips!
 */
PLL_EXPORT pll_split_t * pll_utree_split_create(pll_utree_t * tree,
                                                unsigned int n_tips,
                                                unsigned int * _n_splits)
{
  unsigned int i;
  unsigned int n_splits, split_len, split_size;
  pll_split_t * split_list;
  pll_split_t splits;

  /* as many non-trivial splits as inner branches */
  n_splits   = n_tips - 3;
  split_size = (sizeof(pll_split_base_t) * 8);
  split_len  = (n_tips / split_size) + (n_tips % (sizeof(pll_split_base_t) * 8) > 0);

  split_list = (pll_split_t *) malloc(n_splits * sizeof(pll_split_t));

  if (!split_list)
  {
    return NULL;
  }

  splits = (pll_split_t) calloc(n_splits * split_len, sizeof(pll_split_base_t));

  if (!splits)
  {
    free (split_list);
    return NULL;
  }

  for (i=0; i<n_splits; ++i)
    split_list[i] = splits + i*split_len;

  struct cb_split_params split_data;
  split_data.split_len   = split_len;
  split_data.split_size  = split_size;
  split_data.splits      = split_list;
  split_data.tip_count   = n_tips;

  /* traverse for computing the scripts */
  pll_utree_traverse_apply(tree,
                           0,
                           &cb_get_splits,
                           &split_data);

  /* normalize the splits such that first position is set */
  for (i=0; i<n_splits;++i)
    normalize_split(split_list[i], n_tips);

  /* sort splits and keep first position pointing to the allocated array */
  pll_split_t aux = split_list[0];

  qsort(split_list, n_splits, sizeof(pll_split_t), _cmp_splits);

#ifdef ULTRADEBUG
  /* check */
  for (i=1; i<n_splits; ++i)
    if (compare_splits(split_list[i], split_list[i-1], split_len) <= 0)
    {
      printf("Error %d %d\n", i, _cmp_splits(&split_list[i], &split_list[i-1]));
      pll_utree_show_split(split_list[i], n_tips);
      pll_utree_show_split(split_list[i-1], n_tips);
    }
#endif

  for (i=0; split_list[i] != aux && i < n_splits; ++i);
  assert(i < n_splits);
  if (i>0)
  {
    /* swap */
    void * aux_mem = malloc(sizeof(pll_split_base_t) * split_len);
    memcpy(aux_mem, aux,           sizeof(pll_split_base_t) * split_len);
    memcpy(aux, split_list[0],     sizeof(pll_split_base_t) * split_len);
    memcpy(split_list[0], aux_mem, sizeof(pll_split_base_t) * split_len);
    free(aux_mem);
    split_list[i] = split_list[0];
    split_list[0] = aux;
  }

  if (_n_splits) *_n_splits = n_splits;
  return split_list;
}

PLL_EXPORT void pll_utree_split_destroy(pll_split_t * split_list)
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
 * at positions given by p-matrix index
 */
static int cb_get_splits(pll_utree_t * node, void *data)
{
  struct cb_split_params * split_data = (struct cb_split_params *) data;
  pll_split_t current_split;

  unsigned int tip_count     = split_data->tip_count;
  unsigned int split_size    = split_data->split_size;
  unsigned int split_len     = split_data->split_len;
  unsigned int my_index      = node->pmatrix_index - tip_count,
           child_index;
  unsigned int tip_id, split_id;

  if (!(pll_utree_is_tip(node) || node->pmatrix_index < tip_count))
  {
    /**
     * If the assertion below fails, probably the p-matrix assignments is not
     * suitable for this function.
     */
    assert(my_index < (tip_count - 3));

    current_split = split_data->splits[my_index];

    /* add left */
    if (!pll_utree_is_tip(node->next->back))
    {
      child_index   = node->next->pmatrix_index - tip_count;
      assert(node->next->pmatrix_index >= tip_count && child_index < (tip_count - 3));
      clone_split(current_split, split_data->splits[child_index], split_len);
    }
    else
    {
      tip_id     = node->next->pmatrix_index;
      assert(tip_id < tip_count);
      split_id   = tip_id / split_size;
      tip_id    %= split_size;
      current_split[split_id] = (1 << tip_id);
    }

    /* add right */
    if (!pll_utree_is_tip(node->next->next->back))
    {
      child_index   = node->next->next->pmatrix_index - tip_count;
      assert(node->next->next->pmatrix_index >= tip_count && child_index < (tip_count - 3));
      merge_split(current_split, split_data->splits[child_index], split_len);
    }
    else
    {
      tip_id     = node->next->next->pmatrix_index;
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
      return s1[i] - s2[i];
  }
  return 0;
}

/*
 * Precondition: splits *must* be different.
 */
static int _cmp_splits (const void * a, const void * b)
{
  pll_split_t * s1 = (pll_split_t *) a;
  pll_split_t * s2 = (pll_split_t *) b;
  unsigned int limit = 10000; /* max_taxa = split_size * 10^4 */
  int i = 0;
  for(;((*s1)[i]==(*s2)[i]) && limit;--limit,++i);
  assert(limit);
  return((*s1)[i]-(*s2)[i]);
}

static void normalize_split(pll_split_t split, unsigned int n_tips)
{
  unsigned int split_size = sizeof(pll_split_base_t) * 8;
  unsigned int split_offset = n_tips % split_size;
  unsigned int split_len  = n_tips / split_size + (split_offset>0);
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
