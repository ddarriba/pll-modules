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

kernel_weight_t suggest_weights(const msa_t &msa) {
  kernel_weight_t kw;
  double taxa = msa.count();
  double sites = msa.length();
  double states = msa.states();

  kw[test_kernel_t::partial] =
      0.007156765 * taxa + -2.719444e-05 * sites + 1.328822 + states + 37.70555;
  kw[test_kernel_t::likelihood] = 0.0001289232 * sites + -0.004789004 * taxa +
                                  -0.9559354 * states + 35.20843;
  kw[test_kernel_t::derivative] = 2.240574e-06 * sites + -0.003947451 * taxa +
                                  -0.6615589 * states + 25.72586;
  kw[test_kernel_t::pmatrix] = -0.0001090363 * sites + 0.001506259 * sites +
                               0.2866474 * states + 1.448459;

  for (auto &kv : kw) {
    kv.second = kv.second < 0.0 ? 0.0 : kv.second;
  }
  return kw;
}

kernel_weight_t suggest_weights_2(const msa_t &msa) {
  kernel_weight_t kw;
  double taxa = msa.count();
  double sites = msa.length();
  double states = msa.states();

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

attributes_t select_kernel(const pll_partition_t *pll_partition,
                           const pll_msa_t *pll_msa, const kernel_weight_t &kw,
                           bool fast) {
  msa_t msa{pll_msa};
  model_t model{pll_partition};
  return select_kernel(model, msa, kw, fast);
}

attributes_time_t select_kernel_verbose(const model_t &model, const msa_t &msa,
                                        const kernel_weight_t &kw, bool fast) {
  return fast ? select_kernel_fast_verbose(model, msa, kw)
              : select_kernel_slow_verbose(model, msa, kw);
}

/*
 * Select the attributes in pairs, so that we can eliminate some of the
 * searching. The order is
 *   1. SIMD
 *   2. (Site repeats, tip pattern, none)
 *   3. Rate scalers
 * for the others, we assume a default set of parameters that are likely too be
 * correct.
 */
attributes_time_t select_kernel_fast_verbose(const model_t &model,
                                             const msa_t &msa,
                                             const kernel_weight_t &kw) {
  attributes_time_t times;
  msa_compressed_t cmsa(msa);

  attributes_t attribs(false, true, true, test_cpu_t::avx);

  for (int i = test_cpu_t::avx2; i >= test_cpu_t::none; --i) {
    attribs.simd = static_cast<test_cpu_t>(i);
    test_case_t tc(attribs);
    times[attribs] = weight_kernel_times(kw, tc.benchmark(msa, model));
    std::cout << "timing " << attribs << ":" << times[attribs].count()
              << std::endl;
  }

  attribs = best_attrib_time(times);

  for (size_t i = 0; i < 2; ++i) {
    attribs.site_repeats = false;
    attribs.pattern_tip = static_cast<bool>(i);
    test_case_t tc(attribs);
    times[attribs] = weight_kernel_times(kw, tc.benchmark(msa, model));
    std::cout << "timing " << attribs << ":" << times[attribs].count()
              << std::endl;
  }

  attribs = best_attrib_time(times);

  for (size_t i = 0; i < 1; ++i) {
    attribs.rate_scalers = static_cast<bool>(i);
    test_case_t tc(attribs);
    times[attribs] = weight_kernel_times(kw, tc.benchmark(msa, model));
    std::cout << "timing " << attribs << ":" << times[attribs].count()
              << std::endl;
  }

  return times;
}

attributes_time_t select_kernel_slow_verbose(const model_t &model,
                                             const msa_t &msa,
                                             const kernel_weight_t &kw) {
  attributes_time_t times;
  msa_compressed_t cmsa(msa);
  for (uint8_t bit_attribs = 0; bit_attribs < 0x4; ++bit_attribs) {
    for (uint8_t simd = test_cpu_t::none; simd <= test_cpu_t::avx2; ++simd) {
      if (static_cast<bool>(bit_attribs & (1 << 0)) &&
          static_cast<bool>(bit_attribs & (1 << 1))) {
        continue;
      }

      attributes_t attribs(static_cast<bool>(bit_attribs & (1 << 0)),
                           static_cast<bool>(bit_attribs & (1 << 1)),
                           static_cast<bool>(0), static_cast<test_cpu_t>(simd));
      test_case_t tc(attribs);
      times[attribs] = weight_kernel_times(kw, tc.benchmark(msa, model));
    }
  }
  return times;
}

attributes_t select_kernel(const model_t &model, const msa_t &msa,
                           const kernel_weight_t &kw, bool fast) {
  auto times = select_kernel_verbose(model, msa, kw, fast);
  return best_attrib_time(times);
}

#ifdef __linux__
std::string build_path(size_t cpu_number) {
  std::string s;
  s += "/sys/devices/system/cpu/cpu";
  int leading_zeros = static_cast<int>(std::floor(std::log10(cpu_number)));
  for (int i = 0; i < leading_zeros; ++i) {
    s += "0";
  }
  s += std::to_string(cpu_number);
  s += "/topology/core_id";
  return s;
}

size_t get_core_id(size_t cpu_number) {
  std::string filename = build_path(cpu_number);
  std::ifstream f(filename.c_str());
  char buf[32];
  f.getline(buf, 32);
  return std::stoul(buf);
}

size_t physical_cpu_count() {
  size_t cpu_max = 0;
  size_t n_cpu = static_cast<size_t>(sysconf(_SC_NPROCESSORS_ONLN));
  for (size_t i = 0; i < n_cpu; ++i) {
    size_t core_id = get_core_id(i);
    cpu_max = core_id > cpu_max ? core_id : cpu_max;
  }
  return cpu_max + 1;
}

#elif _WIN32
size_t physical_cpu_count() { return 0; }
#endif
} // namespace dks
