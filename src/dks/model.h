#ifndef DKS_MODEL_H_
#define DKS_MODEL_H_
#include "tree.h"
#include <memory>
#include <pll.h>
#include <vector>

namespace dks {
typedef std::vector<std::vector<char>> msa_t;
typedef std::vector<unsigned int> msa_weight_t;
class model_t {
public:
  model_t(const msa_t &msa, unsigned int states) : model_t{msa, states, 0} {};

  model_t(size_t tip_count, unsigned int states)
      : model_t{tip_count, states, 0} {};

  model_t(const msa_t &msa, unsigned int states, uint64_t seed)
      : _tree{msa.size(), seed}, _states{states}, _rate_categories{1},
        _prop_invar{0.0}, _subst_params((_states - 1) * (_states - 2), 1.0),
        _frequencies(_states, 1.0 / _states){};

  model_t(size_t tip_count, unsigned int states, uint64_t seed)
      : _tree{tip_count, seed}, _states{states}, _rate_categories{1},
        _prop_invar{0.0}, _subst_params{6, 1.0}, _frequencies{4, .25} {};

  model_t(const pll_partition_t *pll_partition)
      : model_t(pll_partition->tips, pll_partition->states, 0,
                pll_partition->rate_cats, pll_partition->subst_params,
                pll_partition->frequencies, *(pll_partition->prop_invar)){};

  model_t(size_t tip_count, unsigned int states, size_t model_index,
          unsigned int rate_categories, double **subst_params,
          double **frequencies, double pinv)
      : _tree{tree_t(tip_count)}, _states{states},
        _rate_categories{rate_categories}, _prop_invar{pinv},
        _subst_params{subst_params[model_index],
                      subst_params[model_index] +
                          (_states - 1) * (_states - 2)},
        _frequencies{frequencies[model_index],
                     frequencies[model_index] + _states} {};

  model_t(size_t tip_count, unsigned int states, unsigned int rate_categories,
          unsigned int seed)
      : _tree{tip_count, seed}, _states{states},
        _rate_categories{rate_categories}, _prop_invar{0.0},
        _subst_params((_states - 1) * (_states - 2) / 2, 1.0),
        _frequencies(_states, 1.0 / _states) {}

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
  unsigned int _states;
  unsigned int _rate_categories;
  double _prop_invar;
  std::vector<double> _subst_params;
  std::vector<double> _frequencies;
};
} // namespace dks
#endif
