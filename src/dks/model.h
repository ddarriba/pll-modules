#pragma once
#include "tree.h"
#include <memory>
#include <pll.h>
#include <vector>

namespace dks {
class model_t {
public:
  model_t(const msa_t &msa) : model_t{msa, 0} {};

  model_t(size_t tip_count) : model_t{tip_count, 0} {};

  model_t(const msa_t &msa, uint64_t seed)
      : _tree{msa.count(), seed}, _states{msa.states()},
        _subst_params((_states - 1) * (_states - 2), 1.0),
        _frequencies(_states, 1.0 / _states){};

  model_t(size_t tip_count, uint64_t seed)
      : _tree{tip_count, seed}, _subst_params{6, 1.0}, _frequencies{4, .25} {};

  model_t(const pll_partition_t *pll_partition)
      : model_t(pll_partition->tips, pll_partition->states, 0,
                pll_partition->subst_params, pll_partition->frequencies){};

  model_t(size_t tip_count, size_t states, size_t model_index,
          double **subst_params, double **frequencies)
      : _tree{tree_t(tip_count)}, _states{states},
        _subst_params{subst_params[model_index],
                      subst_params[model_index] +
                          (_states - 1) * (_states - 2)},
        _frequencies{frequencies[model_index],
                     frequencies[model_index] + _states} {};

  unsigned int submodels() const;
  unsigned int rate_categories() const;
  uint64_t states() const;
  const double *subst_params_raw() const;
  const std::vector<double> &subst_params() const;
  const double *frequencies_raw() const;
  const std::vector<double> &frequencies() const;
  const tree_t &tree() const;

  void reset_tree();

  std::vector<pll_operation_t> make_operations() const;

private:
  tree_t _tree;
  size_t _states;
  std::vector<double> _subst_params;
  std::vector<double> _frequencies;
};
} // namespace dks
