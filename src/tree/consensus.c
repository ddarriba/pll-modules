#include "pll_tree.h"

static void reverse_split(pll_split_t split, unsigned int n_tips);
static int is_subsplit(pll_split_t child,
                       pll_split_t parent,
                       unsigned int split_len);
static pll_utree_t * create_consensus_node(pll_utree_t * parent,
                                           pll_split_t split,
                                           unsigned int split_len,
                                           unsigned int n_tips);
static pll_utree_t * find_splitnode_recurse(pll_split_t split,
                                            pll_utree_t * root,
                                            unsigned int split_len);
static void reset_template_indices(pll_utree_t * node,
                                   unsigned int tip_count);

static pll_split_t clone_split(const pll_split_t from,
                               unsigned int split_len)
{
  pll_split_t to = (pll_split_t) calloc(split_len, sizeof(pll_split_base_t));
  memcpy(to, from, sizeof(pll_split_base_t) * split_len);

  return to;
}

typedef struct consensus_data
{
  pll_split_t split;
  int degree;
  unsigned int n_tips;
  unsigned int split_len;
  unsigned int bit_count;
} consensus_data_t;

static void recursive_assign_indices(pll_utree_t * node,
                                    unsigned int * inner_clv_index,
                                    int * inner_scaler_index,
                                    unsigned int * inner_node_index)
{
  if (!node->next)
  {
    node->clv_index = node->node_index;
    node->pmatrix_index = node->node_index;
    node->scaler_index = PLL_SCALE_BUFFER_NONE;
    return;
  }

  recursive_assign_indices(node->next->back,
                           inner_clv_index,
                           inner_scaler_index,
                           inner_node_index);

  recursive_assign_indices(node->next->next->back,
                           inner_clv_index,
                           inner_scaler_index,
                           inner_node_index);

  node->node_index = *inner_node_index;
  node->next->node_index = *inner_node_index + 1;
  node->next->next->node_index = *inner_node_index + 2;

  node->clv_index = *inner_clv_index;
  node->next->clv_index = *inner_clv_index;
  node->next->next->clv_index = *inner_clv_index;

  node->pmatrix_index = *inner_clv_index;
  node->next->pmatrix_index = node->next->back->pmatrix_index;
  node->next->next->pmatrix_index = node->next->next->back->pmatrix_index;

  node->scaler_index = *inner_scaler_index;
  node->next->scaler_index = *inner_scaler_index;
  node->next->next->scaler_index = *inner_scaler_index;

  *inner_clv_index = *inner_clv_index + 1;
  *inner_scaler_index = *inner_scaler_index + 1;
  *inner_node_index = *inner_node_index + 3;
}

static void reset_template_indices(pll_utree_t * node,
                                   unsigned int tip_count)
{
    unsigned int inner_clv_index = tip_count;
  unsigned int inner_node_index = tip_count;
  int inner_scaler_index = 0;

  recursive_assign_indices(node->back,
                           &inner_clv_index,
                           &inner_scaler_index,
                           &inner_node_index);

  recursive_assign_indices(node->next->back,
                           &inner_clv_index,
                           &inner_scaler_index,
                           &inner_node_index);

  recursive_assign_indices(node->next->next->back,
                           &inner_clv_index,
                           &inner_scaler_index,
                           &inner_node_index);

  node->node_index = inner_node_index;
  node->next->node_index = inner_node_index + 1;
  node->next->next->node_index = inner_node_index + 2;

  node->clv_index = inner_clv_index;
  node->next->clv_index = inner_clv_index;
  node->next->next->clv_index = inner_clv_index;

  node->scaler_index = inner_scaler_index;
  node->next->scaler_index = inner_scaler_index;
  node->next->next->scaler_index = inner_scaler_index;

  node->pmatrix_index = node->back->pmatrix_index;
  node->next->pmatrix_index = node->next->back->pmatrix_index;
  node->next->next->pmatrix_index = node->next->next->back->pmatrix_index;
}

static unsigned int setbit_count(pll_split_t split,
                        unsigned int split_len)
{
  unsigned int setb = 0;
  unsigned int i;

  for (i=0; i<split_len; ++i)
  {
    setb += (unsigned int) __builtin_popcount(split[i]);
  }
  return setb;
}

static unsigned int * get_tip_ids(pll_split_t split,
                                  unsigned int split_len,
                                  unsigned int *count)
{
  unsigned int n_bits = setbit_count(split, split_len);
  unsigned int * ids = (unsigned int *) malloc(n_bits * sizeof(int));
  unsigned int i, base_id, id, pos = 0;
  unsigned int taxa_per_split = 8 * sizeof(split_len);

  for (i=0; i<split_len; ++i)
  {
    base_id = i * taxa_per_split;
    unsigned int s = split[i];
    while(s)
    {
      id = (unsigned int) __builtin_ctz(s);
      if (id < taxa_per_split)
      {
        ids[pos++] = base_id + id;
        s &= ~(1 << id);
      }
    }
  }

  *count = n_bits;
  return ids;
}

static void build_tips_recurse(pll_utree_t * tree)
{
  consensus_data_t * data = (consensus_data_t *) tree->data;
  pll_utree_t * next_root;
  int i;

  next_root = tree;
  for (i=1; i<data->degree; ++i)
  {
    next_root = next_root->next;
    build_tips_recurse(next_root->back);
  }

  //TODO: Temporary only for strict binary
  assert(data->degree == 1 || data->degree == 3);

  if (data->degree == 1)
  {
    /* create tips */
    unsigned int n_bits;
    unsigned int * tip_ids = get_tip_ids(data->split,
                                        data->split_len,
                                        &n_bits);

    //TODO: Temporary only for strict binary
    assert(n_bits == 1);
    if (n_bits == 1)
    {
      tree->node_index = tip_ids[0];
      tree->next = NULL;
    }
    free(tip_ids);
  }
}

static pll_utree_t * find_splitnode_recurse(pll_split_t split,
                                            pll_utree_t * root,
                                            unsigned int split_len)
{
  pll_utree_t * next_root, * ret_node;
  consensus_data_t * data = (consensus_data_t *) root->data;
  int i;

  if (is_subsplit(split, data->split, split_len))
  {
    /* check children */
    next_root = root;
    for (i=1; i<data->degree; ++i)
    {
      next_root = next_root->next;
      if ((ret_node = find_splitnode_recurse(split,
                                 next_root->back,
                                 split_len)) != NULL)
      {
        return ret_node;
      }
    }
    return root;
  }
  else
  {
    return NULL;
  }
}

PLL_EXPORT pll_utree_t * pllmod_utree_from_splits(pll_split_t * splits,
                                                  unsigned int count,
                                                  unsigned int n_tips)
{
  unsigned int split_size = sizeof(pll_split_base_t) * 8;
  unsigned int split_offset = n_tips % split_size;
  unsigned int split_len  = n_tips / split_size + (split_offset>0);
  unsigned int i;
  pll_utree_t * tree, * next_parent;
  pll_split_t rootsplit1, rootsplit2, next_split;

  pll_split_t * all_splits = (pll_split_t *) malloc((count + n_tips) * sizeof(pll_split_t));
  memcpy(all_splits, splits, count * sizeof(pll_split_t));
  for (i=0; i<n_tips; ++i)
  {
    all_splits[count+i] = (pll_split_t) calloc(split_len, sizeof(pll_split_base_t));
    {
      unsigned int split_id   = i / split_size;
      all_splits[count+i][split_id] = (1 << (i % split_size));
    }
  }

  /* create initial tree */
  rootsplit1 = clone_split(splits[0], split_len);
  rootsplit2 = clone_split(splits[0], split_len);
  reverse_split(rootsplit2, n_tips);

  tree = create_consensus_node(NULL, rootsplit1, split_len, n_tips);
  tree->back = create_consensus_node(NULL, rootsplit2, split_len, n_tips);
  tree->back->back = tree;

  /* add splits individually */
  for (i=1; i<(count+n_tips); ++i)
  {
    next_split = clone_split(all_splits[i], split_len);

    /* select branch */
    if (is_subsplit(next_split, rootsplit1, split_len))
      next_parent = tree;
    else if (is_subsplit(next_split, rootsplit2, split_len))
      next_parent = tree->back;
    else
    {
      reverse_split(next_split, n_tips);
      if (is_subsplit(next_split, rootsplit1, split_len))
        next_parent = tree;
      else if (is_subsplit(next_split, rootsplit2, split_len))
        next_parent = tree->back;
      else
      {
        //TODO SET ERROR: Incompatible splits!
        assert(0);
        return NULL;
      }
    }

    next_parent = find_splitnode_recurse(next_split, next_parent, split_len);
    assert(next_parent);

    /* create new node for the split*/
    create_consensus_node(next_parent, next_split, split_len, n_tips);
  }

  build_tips_recurse(tree);
  build_tips_recurse(tree->back);

  reset_template_indices(tree, n_tips);

  return tree;
}












static void connect_consensus_node(pll_utree_t * parent,
                                   pll_utree_t * child,
                                   int auto_rearrange)
{
  consensus_data_t * data_p, * data_c, * data_aux;
  pll_utree_t * new_node, * aux_node, * aux_node2;
  data_p = (consensus_data_t *) parent->data;
  data_c = (consensus_data_t *) child->data;

  assert(is_subsplit(data_c->split, data_p->split, data_p->split_len));
  if (child->back)
  {
    new_node = child->back;

    /* disconnect from parent */
    aux_node = new_node;
    consensus_data_t * data_cp = (consensus_data_t *) new_node->data;
    assert(data_cp->degree > 1);
    while(aux_node->next != child->back) aux_node = aux_node->next;

    --data_cp->degree;
    aux_node->next = new_node->next;
    new_node->next = 0;
  }
  else
  {
    new_node = (pll_utree_t *) malloc (sizeof(pll_utree_t));
  }

  if (auto_rearrange)
  {
    aux_node = parent->next;
    while (aux_node != parent)
    {
      aux_node2 = aux_node->next;
      data_aux = (consensus_data_t *) aux_node->back->data;
      if (is_subsplit(data_aux->split, data_c->split, data_c->split_len))
      {
        connect_consensus_node(child, aux_node->back, 0);
      }
      aux_node = aux_node2;
    }
  }

  /* increase degree */
  ++data_p->degree;

  /* connect new node */
  new_node->data = data_p;
  new_node->next = parent->next;
  parent->next = new_node;

  /* connect child */
  new_node->back = child;
  child->back = new_node;

}

static pll_utree_t * create_consensus_node(pll_utree_t * parent,
                                           pll_split_t split,
                                           unsigned int split_len,
                                           unsigned int n_tips)
{
  pll_utree_t * new_node = (pll_utree_t *) malloc (sizeof(pll_utree_t));
  consensus_data_t * data = (consensus_data_t *)
                                    malloc (sizeof(consensus_data_t));
  data->split     = split;
  data->degree    = 1;
  data->n_tips    = n_tips;
  data->split_len = split_len;

  new_node->data  = data;
  new_node->label = NULL;
  new_node->back = NULL;

  /* self link */
  new_node->next = new_node;

  if (parent)
  {
    connect_consensus_node(parent, new_node, 1);
  }

  return new_node;
}

static int is_subsplit(pll_split_t child,
                       pll_split_t parent,
                       unsigned int split_len)
{
  unsigned int i;
  for (i=0; i<split_len; ++i)
  {
    if ((child[i] & parent[i]) != child[i])
      return 0;
  }
  return 1;
}

static void reverse_split(pll_split_t split, unsigned int n_tips)
{
  unsigned int split_size = sizeof(pll_split_base_t) * 8;
  unsigned int split_offset = n_tips % split_size;
  unsigned int split_len  = n_tips / split_size + (split_offset>0);
  unsigned int i;

  for (i=0; i<split_len; ++i)
    split[i] = ~split[i];

  unsigned int mask = (1<<split_offset) - 1;
  split[split_len - 1] &= mask;
}
