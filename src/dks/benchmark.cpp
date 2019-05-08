#include "benchmark.h"
#include "partition.h"
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <utility>

namespace dks {

msa_t convert_pll_msa_t(const pll_msa_t *pll_msa) {
  msa_t msa;
  msa.reserve(pll_msa->count);
  for (int i = 0; i < pll_msa->count; i++) {
    msa.emplace_back(pll_msa->length);
    for (int j = 0; j < pll_msa->length; j++) {
      msa[i][j] = pll_msa->sequence[i][j];
    }
  }
  return msa;
}

inline benchmark_time_t weight_kernel_times(kernel_weight_t kw,
                                            benchmark_result_t bmr) {
  return kw[test_kernel_t::partial] * bmr[test_kernel_t::partial] +
         kw[test_kernel_t::likelihood] * bmr[test_kernel_t::likelihood] +
         kw[test_kernel_t::derivative] * bmr[test_kernel_t::derivative] +
         kw[test_kernel_t::pmatrix] * bmr[test_kernel_t::pmatrix];
}

inline attributes_t best_attrib_time(const attributes_time_t &at) {
  return std::max_element(at.begin(), at.end(),
                          [](const attributes_time_t::value_type &a,
                             const attributes_time_t::value_type &b) {
                            return a.second > b.second;
                          })
      ->first;
}

kernel_weight_t suggest_weights(double sites, double states, double taxa) {
  kernel_weight_t kw;

  kw[test_kernel_t::partial] =
      0.4866 * sites + 437.1470 * states + 6.5094 * taxa - 3557.8645;
  kw[test_kernel_t::likelihood] =
      0.327 * sites + 28.952 * states + 1.147 * taxa + -43.042;
  kw[test_kernel_t::derivative] =
      0.2174 * sites + 26.1509 * states + 0.7898 * taxa + -108.6298;
  kw[test_kernel_t::pmatrix] =
      3.221e-03 * sites + 2.672e+01 * states + 7.741e-01 * taxa + -2.195e+02;

  for (auto &kv : kw) {
    kv.second = kv.second < 0.0 ? 0.1 : kv.second;
  }
  return kw;
}

kernel_weight_t suggest_weights(const msa_t &msa, unsigned int states) {
  return suggest_weights(msa.size(), msa[0].size(), states);
}

kernel_weight_t suggest_weights(const pllmod_msa_stats_t *msa, int sites,
                                int taxa) {
  return suggest_weights(sites, taxa, msa->states);
}

attributes_time_t
select_kernel_verbose(const model_t &model,
                      const std::vector<std::vector<double>> &clvs,
                      const msa_weight_t &weights, const kernel_weight_t &kw,
                      attributes_generator_t att_gen) {
  attributes_time_t times;
  for (attributes_t attribs = att_gen.next(); attribs != att_gen.end();
       attribs = att_gen.next()) {

    test_case_t tc(attribs);
    times[attribs] =
        weight_kernel_times(kw, tc.benchmark(clvs, weights, model));
  }
  return times;
}

attributes_time_t select_kernel_verbose(const model_t &model, const msa_t &msa,
                                        const msa_weight_t &weights,
                                        const pll_state_t *charmap,
                                        const kernel_weight_t &kw,
                                        attributes_generator_t att_gen) {
  attributes_time_t times;
  for (attributes_t attribs = att_gen.next(); attribs != att_gen.end();
       attribs = att_gen.next()) {

    test_case_t tc(attribs, charmap);
    times[attribs] = weight_kernel_times(kw, tc.benchmark(msa, weights, model));
  }
  return times;
}

unsigned int select_kernel_auto(const pll_partition_t *pll_partition,
                                const pll_msa_t *pll_msa,
                                const pll_state_t *charmap,
                                const kernel_weight_t &kw,
                                attributes_generator_t gen) {
  auto msa = convert_pll_msa_t(pll_msa);
  model_t model{pll_partition};
  msa_weight_t weights(pll_partition->pattern_weights,
                       pll_partition->pattern_weights + pll_msa->length);
  return best_attrib_time(
             select_kernel_verbose(model, msa, weights, charmap, kw, gen))
      .pll_attributes();
}

unsigned int select_kernel_auto(const pll_partition_t *pll_partition,
                                const pll_msa_t *pll_msa,
                                const pll_state_t *charmap) {
  auto kw =
      suggest_weights(pll_msa->length, pll_msa->count, pll_partition->states);
  attributes_generator_t gen;
  return select_kernel_auto(pll_partition, pll_msa, charmap, kw, gen);
}

unsigned int select_kernel_auto(const pll_partition_t *pll_partition,
                                const pll_msa_t *pll_msa,
                                const pll_state_t *charmap,
                                attributes_generator_t gen) {
  auto kw =
      suggest_weights(pll_msa->length, pll_msa->count, pll_partition->states);
  return select_kernel_auto(pll_partition, pll_msa, charmap, kw, gen);
}

unsigned int select_kernel_auto(const msa_t &msa, const msa_weight_t &weights,
                                const pll_state_t *charmap, unsigned int states,
                                unsigned int rate_cats,
                                attributes_generator_t gen) {
  auto kw = suggest_weights(weights.size(), msa.size(), states);
  model_t model{msa, states, rate_cats};
  auto result = select_kernel_verbose(model, msa, weights, charmap, kw, gen);
  return best_attrib_time(result).pll_attributes();
}

unsigned int select_kernel_auto(const msa_t &msa, const msa_weight_t &weights,
                                const pll_state_t *charmap, unsigned int states,
                                unsigned int rate_cats) {
  attributes_generator_t gen;
  return select_kernel_auto(msa, weights, charmap, states, rate_cats, gen);
}

unsigned int select_kernel_auto(const std::vector<std::vector<double>> &clvs,
                                const msa_weight_t &weights,
                                unsigned int states, unsigned int rate_cats,
                                attributes_generator_t gen) {
  auto kw = suggest_weights(weights.size(), clvs.size(), states);
  model_t model{clvs.size(), states, rate_cats};
  auto result = select_kernel_verbose(model, clvs, weights, kw, gen);
  return best_attrib_time(result).pll_attributes();
}

unsigned int select_kernel_auto(const std::vector<std::vector<double>> &clvs,
                                const msa_weight_t &weights,
                                unsigned int states, unsigned int rate_cats) {
  attributes_generator_t gen;
  return select_kernel_auto(clvs, weights, states, rate_cats, gen);
}
} // namespace dks
