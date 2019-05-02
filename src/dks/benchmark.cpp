#include "benchmark.h"
#include "msa.h"
#include "partition.h"
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <utility>

#ifdef __linux__
#include <fstream>
#include <unistd.h>
#elif _WIN32
// nothing
#endif

namespace dks {

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

kernel_weight_t suggest_weights(const msa_t &msa) {
  return suggest_weights(msa.count(), msa.length(), msa.states());
}

kernel_weight_t suggest_weights(const pllmod_msa_stats_t *msa, int sites,
                                int taxa) {
  return suggest_weights(sites, taxa, msa->states);
}

attributes_t select_kernel(const pll_partition_t *pll_partition,
                           const pll_msa_t *pll_msa,
                           const kernel_weight_t &kw) {
  msa_t msa{pll_msa};
  model_t model{pll_partition};
  return select_kernel(model, msa, kw);
}

attributes_t select_kernel_auto(const pll_partition_t *pll_partition,
                                const pll_msa_t *pll_msa) {
  auto kw = suggest_weights(pll_msa);
  return select_kernel(pll_partition, pll_msa, kw);
}

attributes_time_t select_kernel_verbose(const model_t &model, const msa_t &msa,
                                        const kernel_weight_t &kw,
                                        attributes_generator_t att_gen) {
  attributes_time_t times;
  for (attributes_t attribs = att_gen.next(); attribs != att_gen.end();
       attribs = att_gen.next()) {

    test_case_t tc(attribs);
    times[attribs] = weight_kernel_times(kw, tc.benchmark(msa, model));
  }
  return times;
} // namespace dks

attributes_t select_kernel(const model_t &model, const msa_t &msa,
                           const kernel_weight_t &kw,
                           attributes_generator_t gen) {
  auto times = select_kernel_verbose(model, msa, kw, gen);
  return best_attrib_time(times);
}

attributes_t select_kernel(const model_t &model, const msa_t &msa,
                           const kernel_weight_t &kw) {
  attributes_generator_t gen;
  return select_kernel(model, msa, kw, gen);
}
} // namespace dks
