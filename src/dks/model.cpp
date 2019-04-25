#include "model.h"
#include <pll.h>

namespace dks {

unsigned int model_t::submodels() const { return 1; }

unsigned int model_t::rate_categories() const { return 1; }

uint64_t model_t::states() const { return _states; }

const double *model_t::subst_params_raw() const { return _subst_params.data(); }

const std::vector<double> &model_t::subst_params() const {
  return _subst_params;
}

const double *model_t::frequencies_raw() const { return _frequencies.data(); }

const std::vector<double> &model_t::frequencies() const { return _frequencies; }

const tree_t &model_t::tree() const { return _tree; }

std::vector<pll_operation_t> model_t::make_operations() const {
  auto traversal_nodes = _tree.full_traverse();
  std::vector<pll_operation_t> operations(traversal_nodes.size());
  unsigned int operations_count = 0;
  pll_utree_create_operations(traversal_nodes.data(), traversal_nodes.size(),
                              nullptr, nullptr, operations.data(), nullptr,
                              &operations_count);

  operations.resize(operations_count);
  return operations;
}

void model_t::reset_tree() { _tree = tree_t(_tree.tip_count()); }
} // namespace dks
