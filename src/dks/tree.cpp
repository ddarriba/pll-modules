#include "pll.h"
#include "tree.h"
#include <random>

namespace dks {
int full_traverse_cb(pll_unode_t *n) {
  // silence a warning
  do {
    (void)(n);
  } while (0);
  return PLL_SUCCESS;
}

tree_t::tree_t(size_t tip_count, uint64_t random_seed) {
  size_t inner_count = tip_count - 2;

  pll_unode_t *vroot = make_triplet(make_tip(), make_tip(), make_tip());

  // generate a list of nodes
  size_t node_buffer_size = tip_count + 3 * inner_count;
  pll_unode_t **nodes =
      (pll_unode_t **)malloc(sizeof(pll_unode_t *) * node_buffer_size);

  std::minstd_rand random_engine(random_seed);

  for (size_t i = 3; i < tip_count; i++) {
    unsigned int node_count;
    pll_utree_traverse(vroot, PLL_TREE_TRAVERSE_POSTORDER, full_traverse_cb,
                       nodes, &node_count);

    std::uniform_int_distribution<> roller(0, node_count - 1);
    pll_unode_t *insert_node = nodes[roller(random_engine)];
    insert_tip(insert_node);
  }
  _tree = pll_utree_wraptree(vroot, tip_count);
  pll_utree_reset_template_indices(vroot, tip_count);
  update_internal_lists();

  free(nodes);
  if (!pll_utree_check_integrity(_tree)) {
    pll_utree_destroy(_tree, nullptr);
    throw std::runtime_error("Tree validation check failed");
  }
}

tree_t::~tree_t() { pll_utree_destroy(_tree, nullptr); }

pll_unode_t *tree_t::vroot() const { return _tree->vroot; }

size_t tree_t::node_count() const {
  return _tree->tip_count + _tree->inner_count;
}

size_t tree_t::edge_count() const { return _tree->tip_count * 2 - 3; }

size_t tree_t::tip_count() const { return _tree->tip_count; }

void tree_t::insert_tip(pll_unode_t *insert_node) {
  pll_unode_t *new_a = make_node();
  pll_unode_t *new_b = make_node();
  pll_unode_t *new_c = make_node();
  pll_unode_t *new_tip = make_tip();

  pll_unode_t *saved_back = insert_node->back;

  pair_nodes(insert_node, new_a);
  pair_nodes(saved_back, new_b);
  pair_nodes(new_tip, new_c);

  make_circle(new_a, new_b, new_c);
}

void tree_t::pair_nodes(pll_unode_t *a, pll_unode_t *b) {
  a->back = b;
  b->back = a;
}

void tree_t::make_circle(pll_unode_t *a, pll_unode_t *b, pll_unode_t *c) {
  a->next = b;
  b->next = c;
  c->next = a;
}

pll_unode_t *tree_t::make_triplet(pll_unode_t *a, pll_unode_t *b,
                                  pll_unode_t *c) {
  pll_unode_t *inner_a = make_node();
  pll_unode_t *inner_b = make_node();
  pll_unode_t *inner_c = make_node();

  pair_nodes(inner_a, a);
  pair_nodes(inner_b, b);
  pair_nodes(inner_c, c);

  make_circle(inner_a, inner_b, inner_c);
  return inner_a;
}

pll_unode_t *tree_t::make_node() {
  pll_unode_t *a = (pll_unode_t *)malloc(sizeof(pll_unode_t));
  a->label = nullptr;
  a->next = nullptr;
  a->back = nullptr;
  a->length = 0.1;
  return a;
}

pll_unode_t *tree_t::make_tip() { return make_node(); }

const std::vector<double> &tree_t::branch_lengths() const {
  return _branch_lengths;
}

const std::vector<unsigned int> &tree_t::matrix_indices() const {
  return _matrix_indices;
}

std::vector<pll_unode_t *> tree_t::full_traverse() const {
  unsigned int node_number = 0;
  std::vector<pll_unode_t *> trav_nodes(node_count());
  if (!(pll_utree_traverse(_tree->vroot, PLL_TREE_TRAVERSE_POSTORDER,
                           full_traverse_cb, trav_nodes.data(),
                           &node_number))) {
    throw std::runtime_error("Generic failure to traverse the tree");
  }

  return trav_nodes;
}

void tree_t::fill_branch_lengths(const std::vector<pll_unode_t *> &nodes) {
  if (_branch_lengths.size() != nodes.size()) {
    _branch_lengths.resize(nodes.size(), 0.1);
  }
  for (size_t i = 0; i < nodes.size(); i++) {
    _branch_lengths[i] = nodes[i]->length;
  }
}

void tree_t::fill_matrix_indices(const std::vector<pll_unode_t *> &nodes) {
  if (_matrix_indices.size() != nodes.size()) {
    _matrix_indices.resize(nodes.size(), 0);
  }
  for (size_t i = 0; i < nodes.size(); i++) {
    _matrix_indices[i] = nodes[i]->pmatrix_index;
  }
}

void tree_t::update_internal_lists() {
  auto nodes = full_traverse();
  fill_branch_lengths(nodes);
  fill_matrix_indices(nodes);
}

std::vector<pll_operation_t> tree_t::make_operations() const {
  auto traversal_nodes = full_traverse();
  std::vector<pll_operation_t> operations(traversal_nodes.size());
  unsigned int operations_count = 0;
  pll_utree_create_operations(traversal_nodes.data(), traversal_nodes.size(),
                              nullptr, nullptr, operations.data(), nullptr,
                              &operations_count);

  operations.resize(operations_count);
  return operations;
}
} // namespace dks
