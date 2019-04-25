#pragma once
#include "model.h"
#include "msa.h"
#include <memory>
#include <pll.h>
#include <utility>

namespace dks {
typedef std::pair<double, double> derivative_t;
class partition_t {
public:
  partition_t(unsigned int tips, unsigned int clv_buffers, unsigned int states,
              unsigned int sites, unsigned int rate_matrices,
              unsigned int prob_matrices, unsigned int rate_cats,
              unsigned int scale_buffers, unsigned int attributes)
      : _partition{pll_partition_create(tips, clv_buffers, states, sites,
                                        rate_matrices, prob_matrices, rate_cats,
                                        scale_buffers, attributes)} {};
  partition_t(const msa_t &, const model_t &, unsigned int);
  ~partition_t();

  void initialize_tips(const msa_t &);
  void initialize_rates(const model_t &model);
  void set_pattern_weights(const msa_compressed_t &);
  void set_pattern_weights(const msa_t &);
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
