#pragma once
#include "model.h"
#include <memory>
#include <pll.h>
#include <utility>

namespace dks {
typedef std::pair<double, double> derivative_t;
class partition_t {
public:
  partition_t(const msa_t &, const model_t &, const msa_weight_t &,
              const pll_state_t *, unsigned int);
  partition_t(const std::vector<std::vector<double>> &, const model_t &,
              const msa_weight_t &, unsigned int);
  ~partition_t();

  void initialize_tips(const msa_t &, const pll_state_t *);
  void initialize_tips(const std::vector<std::vector<double>> &);

  void initialize_rates(const model_t &model);
  void set_pattern_weights(const msa_weight_t &);
  void update_probability_matrices(const tree_t &tree);

  void update_partials(const std::vector<pll_operation_t> &);
  void update_partials(const tree_t &tree);
  void update_partials(const model_t &model);

  void update_site_repeats(const std::vector<pll_operation_t> &);
  void update_site_repeats(const tree_t &tree);
  void update_site_repeats(const model_t &model);

  void invalidate_prob_matrix();

  void update_sumtable(const tree_t &tree);
  derivative_t compute_derivative(const tree_t &tree, double brlen = 1.0);

  double loglh(const tree_t &tree);
  double loglh(const model_t &model);
  std::vector<double> loglh_persite(const model_t &model, size_t sites);

  bool site_repeats() const;
  unsigned int attributes() const;

  const pll_partition_t *partition_raw() const;

private:
  void alloc_sumtable(unsigned int attribs);

  pll_partition_t *_partition;
  double *_sumtable;
  constexpr static unsigned int _params_indices[] = {0};
  constexpr static double _rate_cats[] = {1.0};
};
} // namespace dks
