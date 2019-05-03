#include "partition.h"
#include <pll.h>

namespace dks {

constexpr unsigned int partition_t::_params_indices[];
constexpr double partition_t::_rate_cats[];

partition_t::partition_t(const msa_t &msa, const model_t &model,
                         const msa_weight_t &weights,
                         const pll_state_t *charmap, unsigned int attributes) {

  unsigned int tip_count = msa.size();
  unsigned int inner_count = tip_count - 2;

  _partition = pll_partition_create(tip_count,               // tips
                                    inner_count,             // clv_buffers
                                    model.states(),          // states
                                    msa[0].size(),           // sites
                                    model.submodels(),       // rate_matrices
                                    2 * tip_count - 3,       // prob_matrices
                                    model.rate_categories(), // rate_cats
                                    inner_count,             // scale_buffes
                                    attributes               // attributes
  );
  initialize_tips(msa, charmap);
  initialize_rates(model);
  set_pattern_weights(weights);
  update_probability_matrices(model.tree());
  alloc_sumtable(attributes);
}

partition_t::partition_t(const std::vector<std::vector<double>> &clvs,
                         const model_t &model, const msa_weight_t &weights,
                         unsigned int attributes) {

  unsigned int tip_count = clvs.size();
  unsigned int inner_count = tip_count - 2;

  _partition = pll_partition_create(tip_count,               // tips
                                    inner_count,             // clv_buffers
                                    model.states(),          // states
                                    clvs[0].size(),          // sites
                                    model.submodels(),       // rate_matrices
                                    2 * tip_count - 3,       // prob_matrices
                                    model.rate_categories(), // rate_cats
                                    inner_count,             // scale_buffes
                                    attributes               // attributes
  );
  initialize_tips(clvs);
  initialize_rates(model);
  set_pattern_weights(weights);
  update_probability_matrices(model.tree());
  alloc_sumtable(attributes);
}

void partition_t::alloc_sumtable(unsigned int attribs) {
  unsigned int align = PLL_ALIGNMENT_CPU;
  switch (attribs & PLL_ATTRIB_ARCH_MASK) {
  case PLL_ATTRIB_ARCH_SSE:
    align = PLL_ALIGNMENT_SSE;
    break;
  case PLL_ATTRIB_ARCH_AVX:
  case PLL_ATTRIB_ARCH_AVX2:
  case PLL_ATTRIB_ARCH_AVX512:
    align = PLL_ALIGNMENT_AVX;
    break;
  default:
    break;
  }
  _sumtable = (double *)pll_aligned_alloc(
      _partition->sites * _partition->rate_cats * _partition->states_padded *
          sizeof(double),
      align);
}

void partition_t::initialize_tips(const msa_t &msa,
                                  const pll_state_t *charmap) {
  for (size_t tip_id = 0; tip_id < msa.size(); tip_id++) {
    pll_set_tip_states(_partition, tip_id, charmap, msa[tip_id].data());
  }
}

void partition_t::initialize_tips(const std::vector<std::vector<double>> &clvs){
  for (size_t tip_id = 0; tip_id < clvs.size(); tip_id++) {
    pll_set_tip_clv(_partition, tip_id, clvs[tip_id].data(), PLL_FALSE);
  }
}

void partition_t::initialize_rates(const model_t &model) {
  for (size_t i = 0; i < model.submodels(); i++) {
    pll_set_subst_params(_partition, i, model.subst_params_raw());
    pll_set_frequencies(_partition, i, model.frequencies_raw());
  }
  pll_set_category_rates(_partition, _rate_cats);
}

void partition_t::update_probability_matrices(const tree_t &tree) {
  pll_update_prob_matrices(
      _partition, _params_indices, tree.matrix_indices().data(),
      tree.branch_lengths().data(), tree.branch_lengths().size());
}

void partition_t::update_partials(const std::vector<pll_operation_t> &ops) {
  pll_update_partials(_partition, ops.data(), ops.size());
}

void partition_t::update_site_repeats(const tree_t &tree) {
  update_site_repeats(tree.make_operations());
}

void partition_t::update_site_repeats(const model_t &model) {
  update_site_repeats(model.make_operations());
}

void partition_t::update_site_repeats(const std::vector<pll_operation_t> &ops) {
  for (size_t i = 0; i < ops.size(); i++) {
    pll_update_repeats(_partition, &ops[i]);
  }
}

void partition_t::set_pattern_weights(const msa_weight_t &weights) {
  pll_set_pattern_weights(_partition, weights.data());
}

void partition_t::update_partials(const tree_t &tree) {
  update_partials(tree.make_operations());
}

void partition_t::update_partials(const model_t &model) {
  update_partials(model.tree());
}

void partition_t::invalidate_prob_matrix() {
  *_partition->eigen_decomp_valid = 0;
}

double partition_t::loglh(const tree_t &tree) {
  pll_unode_t *parent = tree.vroot();
  pll_unode_t *child = parent->back;
  return pll_compute_edge_loglikelihood(
      _partition, parent->clv_index, parent->scaler_index, child->clv_index,
      child->scaler_index, parent->pmatrix_index, _params_indices, nullptr);
}

double partition_t::loglh(const model_t &model) { return loglh(model.tree()); }

std::vector<double> partition_t::loglh_persite(const model_t &model,
                                               size_t sites) {
  pll_unode_t *parent = model.tree().vroot();
  pll_unode_t *child = parent->back;

  std::vector<double> persites(sites, 0.0);

  pll_compute_edge_loglikelihood(_partition, parent->clv_index,
                                 parent->scaler_index, child->clv_index,
                                 child->scaler_index, parent->pmatrix_index,
                                 _params_indices, persites.data());
  return persites;
}

unsigned int partition_t::attributes() const { return _partition->attributes; }

void partition_t::update_sumtable(const tree_t &tree) {
  pll_unode_t *parent = tree.vroot();
  pll_unode_t *child = parent->back;
  if (site_repeats()) {
    update_partials(tree);
  }

  pll_update_sumtable(_partition, parent->clv_index, child->clv_index,
                      parent->scaler_index, child->scaler_index,
                      _params_indices, _sumtable);
}

derivative_t partition_t::compute_derivative(const tree_t &tree, double brlen) {
  pll_unode_t *parent = tree.vroot();
  pll_unode_t *child = parent->back;

  derivative_t d;
  pll_compute_likelihood_derivatives(
      _partition, parent->scaler_index, child->scaler_index, brlen,
      _params_indices, _sumtable, &d.first, &d.second);

  return d;
}

bool partition_t::site_repeats() const {
  return attributes() & PLL_ATTRIB_SITE_REPEATS;
}

partition_t::~partition_t() {
  pll_partition_destroy(_partition);
  pll_aligned_free(_sumtable);
}

const pll_partition_t *partition_t::partition_raw() const { return _partition; }

} // namespace dks
