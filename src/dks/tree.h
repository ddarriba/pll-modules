#pragma once
#include "pll.h"
#include <vector>
#include <memory>
#include <utility>

namespace dks {
class tree_t {
public:
  tree_t(size_t tip_count) : tree_t{tip_count, 0} {};
  tree_t(size_t tip_count, uint64_t random_seed);
  pll_unode_t *vroot() const;
  size_t node_count() const;
  size_t edge_count() const;
  size_t tip_count() const;
  const std::vector<double> &branch_lengths() const;
  const std::vector<unsigned int> &matrix_indices() const;
  std::vector<pll_unode_t *> full_traverse() const;
  std::vector<pll_operation_t> make_operations() const;
  ~tree_t();

private:
  static void insert_tip(pll_unode_t *);
  static void pair_nodes(pll_unode_t *, pll_unode_t *);
  static void make_circle(pll_unode_t *, pll_unode_t *, pll_unode_t *);
  static pll_unode_t *make_triplet(pll_unode_t *, pll_unode_t *, pll_unode_t *);
  static pll_unode_t *make_node();
  static pll_unode_t *make_tip();

  void fill_branch_lengths(const std::vector<pll_unode_t *> &);
  void fill_matrix_indices(const std::vector<pll_unode_t *> &);

  void update_internal_lists();

  pll_utree_t *_tree;
  std::vector<double> _branch_lengths;
  std::vector<unsigned int> _matrix_indices;
};
} // namespace dks
