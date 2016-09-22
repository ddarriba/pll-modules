#include "pll_tree.h"
#include "../pllmod_common.h"

typedef struct consensus_data
{
  pll_split_t split;
  int degree;
  unsigned int tip_count;
  unsigned int split_len;
  unsigned int bit_count;
} consensus_data_t;

static void reverse_split(pll_split_t split, unsigned int tip_count);
static int is_subsplit(pll_split_t child,
                       pll_split_t parent,
                       unsigned int split_len);
static unsigned int setbit_count(pll_split_t split,
                                 unsigned int split_len);
static pll_utree_t * find_splitnode_recurse(pll_split_t split,
                                            pll_utree_t * root,
                                            unsigned int split_len);
static int get_split_id(pll_split_t split,
                        unsigned int split_len);
static pll_utree_t * create_consensus_node(pll_utree_t * parent,
                                           pll_split_t split,
                                           unsigned int split_len,
                                           unsigned int tip_count);
static pll_utree_t * find_splitnode_recurse(pll_split_t split,
                                            pll_utree_t * root,
                                            unsigned int split_len);
static pll_split_t clone_split(const pll_split_t from,
                               unsigned int split_len);
static void recursive_assign_indices(pll_utree_t * node,
                                     unsigned int * inner_clv_index,
                                     int * inner_scaler_index,
                                     unsigned int * inner_node_index);
static void reset_template_indices(pll_utree_t * node,
                                   unsigned int tip_count);
static void build_tips_recurse(pll_utree_t * tree, const char **tip_labels);


static int cb_destroy_data(pll_utree_t * node, void * d)
{
  UNUSED(d);
  consensus_data_t * data = (consensus_data_t *) node->data;
  free(data->split);
  free(data);
  node->data = NULL;
  return PLL_SUCCESS;
}

PLL_EXPORT int pllmod_utree_compatible_splits(pll_split_t s1,
                                              pll_split_t s2,
                                              unsigned int split_len)
{
  unsigned int i;

  /* check conflicts between s1 and s2 */
  for(i = 0; i < split_len; i++)
    if(s1[i] & s2[i])
      break;

  if(i == split_len)
    return 1;

  /* check conflicts between s1 and ~s2 */
  for(i = 0; i < split_len; i++)
    if(s1[i] & ~s2[i])
      break;

  if(i == split_len)
    return 1;

  /* check conflicts between ~s1 and s2 */
  for(i = 0; i < split_len; i++)
    if(~s1[i] & s2[i])
      break;

  if(i == split_len)
    return 1;
  else
    return 0;
}

PLL_EXPORT pll_utree_t * pllmod_utree_from_splits(const pll_split_t * splits,
                                                  unsigned int split_count,
                                                  unsigned int tip_count,
                                                  const char **tip_labels)
{
  unsigned int split_size = sizeof(pll_split_base_t) * 8;
  unsigned int split_offset = tip_count % split_size;
  unsigned int split_len  = tip_count / split_size + (split_offset>0);
  unsigned int i;
  pll_utree_t * tree, * next_parent, * return_tree;
  pll_split_t rootsplit1, rootsplit2, next_split;

  pll_split_t * all_splits = (pll_split_t *) malloc((split_count + tip_count) *
                                                    sizeof(pll_split_t));
  memcpy(all_splits, splits, split_count * sizeof(pll_split_t));
  for (i=0; i<tip_count; ++i)
  {
    all_splits[split_count+i] = (pll_split_t) calloc(split_len,
                                               sizeof(pll_split_base_t));
    {
      unsigned int split_id   = i / split_size;
      all_splits[split_count+i][split_id] = (1 << (i % split_size));
    }
  }

  /* create initial tree with 2 connected nodes of degree 1 */
  rootsplit1 = clone_split(splits[0], split_len);
  rootsplit2 = clone_split(splits[0], split_len);
  reverse_split(rootsplit2, tip_count);

  /* build tree out of the first split */
  tree = create_consensus_node(NULL, rootsplit1, split_len, tip_count);
  tree->back = create_consensus_node(NULL, rootsplit2, split_len, tip_count);
  tree->back->back = tree;

  return_tree = tree;

  /* add splits individually */
  for (i=1; i<(split_count+tip_count); ++i)
  {
    next_split = clone_split(all_splits[i], split_len);

    /* select branch */
    if (is_subsplit(next_split, rootsplit1, split_len))
      next_parent = tree;
    else if (is_subsplit(next_split, rootsplit2, split_len))
      next_parent = tree->back;
    else
    {
      reverse_split(next_split, tip_count);
      if (is_subsplit(next_split, rootsplit1, split_len))
        next_parent = tree;
      else if (is_subsplit(next_split, rootsplit2, split_len))
        next_parent = tree->back;
      else
      {
        pllmod_set_error(PLLMOD_TREE_ERROR_INVALID_SPLIT,
                         "Splits are incompatible");
        return_tree = NULL;
        break;
      }
    }

    /* select node */
    next_parent = find_splitnode_recurse(next_split, next_parent, split_len);
    assert(next_parent);

    /* create new node for the split*/
    create_consensus_node(next_parent, next_split, split_len, tip_count);
  }

  build_tips_recurse(tree, tip_labels);
  build_tips_recurse(tree->back, tip_labels);

  reset_template_indices(tree, tip_count);

  /* clean */
  for (i=0; i<tip_count; ++i)
    free(all_splits[split_count+i]);
  free(all_splits);
  pllmod_utree_traverse_apply(tree,
                              cb_destroy_data, /* pre  */
                              NULL,            /*  in  */
                              NULL,            /* post */
                              NULL);

  /* return_tree == tree if sucess, or null if the algorithm failed */
  return return_tree;
}







/******************************************************************************/
/* static functions */

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

static int get_split_id(pll_split_t split,
                        unsigned int split_len)
{
  unsigned int n_bits = setbit_count(split, split_len);
  int i, base_id, ctz;
  int taxa_per_split = 8 * sizeof(split_len);
  int id = -1;

  if (n_bits != 1)
  {
    pllmod_set_error(PLLMOD_TREE_ERROR_INVALID_SPLIT,
                    "Invalid trivial split");
    return -1;
  }

  for (i=0; i<(int)split_len; ++i)
  {
    if (split[i])
    {
      base_id = i * taxa_per_split;
      ctz = (int)__builtin_ctz(split[i]);
      assert (ctz < taxa_per_split);
      id = base_id + ctz;
      break;
    }
  }

  /* if the assertion below fails, there is an error either in this algorithm,
     or in setbit_count. */
  assert(id != -1);

  return id;
}

static void build_tips_recurse(pll_utree_t * tree, const char ** tip_labels)
{
  consensus_data_t * data = (consensus_data_t *) tree->data;
  pll_utree_t * next_root;
  int i;

  next_root = tree;
  for (i=1; i<data->degree; ++i)
  {
    next_root = next_root->next;
    build_tips_recurse(next_root->back, tip_labels);
  }

  if (data->degree == 1)
  {
    /* create tips */
    int tip_id = get_split_id(data->split,
                              data->split_len);

    assert (tip_id != -1);

    if (tip_labels)
    {
      tree->label = (char *) malloc(strlen(tip_labels[tip_id])+1);
      strcpy(tree->label, tip_labels[tip_id]);
    }
    else
    {
      tree->label = NULL;
    }
    tree->node_index = (unsigned int) tip_id;
    tree->next = NULL;
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
                                           unsigned int tip_count)
{
  pll_utree_t * new_node = (pll_utree_t *) malloc (sizeof(pll_utree_t));
  consensus_data_t * data = (consensus_data_t *)
                                    malloc (sizeof(consensus_data_t));
  data->split     = split;
  data->degree    = 1;
  data->tip_count    = tip_count;
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

static void reverse_split(pll_split_t split, unsigned int tip_count)
{
  unsigned int split_size = sizeof(pll_split_base_t) * 8;
  unsigned int split_offset = tip_count % split_size;
  unsigned int split_len  = tip_count / split_size + (split_offset>0);
  unsigned int i;

  for (i=0; i<split_len; ++i)
    split[i] = ~split[i];

  unsigned int mask = (1<<split_offset) - 1;
  split[split_len - 1] &= mask;
}

static pll_split_t clone_split(const pll_split_t from,
                               unsigned int split_len)
{
  pll_split_t to = (pll_split_t) calloc(split_len, sizeof(pll_split_base_t));
  memcpy(to, from, sizeof(pll_split_base_t) * split_len);

  return to;
}

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
