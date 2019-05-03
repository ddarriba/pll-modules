#include "partition.h"
#include "test_case.h"
#include <chrono>
#include <exception>
#include <ostream>
#include <random>

namespace dks {

test_cpu_t test_cpu_from_attribs(uint32_t attribs) {
  if ((attribs & PLL_ATTRIB_ARCH_MASK) == 0) {
    return dks::test_cpu_t::none;
  }
  return static_cast<test_cpu_t>(32 -
                                 __builtin_clz(attribs & PLL_ATTRIB_ARCH_MASK));
}

int attributes_t::cpu_attrib() const {
  if (simd == dks::test_cpu_t::avx2) {
    return PLL_ATTRIB_ARCH_AVX2;
  }
  if (simd == dks::test_cpu_t::avx) {
    return PLL_ATTRIB_ARCH_AVX;
  }
  if (simd == dks::test_cpu_t::sse) {
    return PLL_ATTRIB_ARCH_SSE;
  }
  if (simd == dks::test_cpu_t::none) {
    return PLL_ATTRIB_ARCH_CPU;
  }
  if (simd == dks::test_cpu_t::avx512) {
    return PLL_ATTRIB_ARCH_AVX512;
  }
  throw std::runtime_error("Unrecognized CPU type");
}

int attributes_t::siterepeat_attrib() const {
  return site_repeats ? PLL_ATTRIB_SITE_REPEATS : 0;
}

int attributes_t::patterntip_attrib() const {
  return pattern_tip ? PLL_ATTRIB_PATTERN_TIP : 0;
}

int attributes_t::pll_attributes() const {
  return cpu_attrib() | siterepeat_attrib() | patterntip_attrib();
}

attributes_t attributes_generator_t::next() {
  while (!valid()) {
    _state++;
    if (_state > _max) {
      return end();
    }
  }
  attributes_t attrib(_state & PLL_ATTRIB_PATTERN_TIP,
                      _state & PLL_ATTRIB_SITE_REPEATS, 0,
                      test_cpu_from_attribs(_state));
  _state++;
  return attrib;
}

attributes_t attributes_generator_t::end() {
  return attributes_t(false, false, false, dks::test_cpu_t::invalid);
}

benchmark_result_t test_case_t::run_benchmarks(partition_t &partition,
                                               const model_t &model) {
  benchmark_result_t br;

  br[test_kernel_t::partial] = benchmark_partials(partition, model);
  br[test_kernel_t::likelihood] = benchmark_likelihood(partition, model);
  br[test_kernel_t::pmatrix] = benchmark_pmatrix(partition, model);
  br[test_kernel_t::derivative] = benchmark_derivative(partition, model);

  return br;
}

benchmark_result_t
test_case_t::benchmark(const std::vector<std::vector<double>> &clvs ,
                       const msa_weight_t &weights, const model_t &model) {
  partition_t partition(clvs, model, weights, _charmap, attributes());

  return run_benchmarks(partition, model);
}

benchmark_result_t test_case_t::benchmark(const msa_t &msa,
                                          const msa_weight_t &weights,
                                          const model_t &model) {
  partition_t partition(msa, model, weights, _charmap, attributes());

  return run_benchmarks(partition, model);
}

benchmark_time_t test_case_t::benchmark_partials(partition_t &partition,
                                                 const model_t &model) {
  for (size_t i = 0; i < _trials; i++) {
    partition.update_partials(model);
  }
  auto t1 = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < _trials; i++) {
    partition.update_partials(model);
  }
  auto t2 = std::chrono::high_resolution_clock::now();
  return (t2 - t1) / _trials;
}

benchmark_time_t test_case_t::benchmark_likelihood(partition_t &partition,
                                                   const model_t &model) {
  partition.update_partials(model);
  for (size_t i = 0; i < _trials; i++) {
    partition.loglh(model);
  }
  auto t1 = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < _trials; i++) {
    partition.loglh(model);
  }
  auto t2 = std::chrono::high_resolution_clock::now();
  return (t2 - t1) / _trials;
}

benchmark_time_t test_case_t::benchmark_pmatrix(partition_t &partition,
                                                const model_t &model) {
  for (size_t i = 0; i < _trials; i++) {
    partition.update_probability_matrices(model.tree());
  }
  auto t1 = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < _trials; i++) {
    partition.update_probability_matrices(model.tree());
  }
  auto t2 = std::chrono::high_resolution_clock::now();
  return (t2 - t1) / _trials;
}

benchmark_time_t test_case_t::benchmark_derivative(partition_t &partition,
                                                   const model_t &model) {
  for (size_t i = 0; i < _trials; i++) {
    partition.update_sumtable(model.tree());
    partition.compute_derivative(model.tree());
  }
  auto t1 = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < _trials; i++) {
    partition.update_sumtable(model.tree());
    partition.compute_derivative(model.tree());
  }
  auto t2 = std::chrono::high_resolution_clock::now();
  return (t2 - t1) / _trials;
}

benchmark_time_t
test_case_t::benchmark_update_site_repeats(partition_t &partition,
                                           const model_t &model) {
  for (size_t i = 0; i < _trials; i++) {
    partition.update_site_repeats(model);
  }
  auto t1 = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < _trials; i++) {
    partition.update_site_repeats(model);
  }
  auto t2 = std::chrono::high_resolution_clock::now();
  return (t2 - t1) / _trials;
}

attributes_t test_case_t::attributes_struct() const {
  return attributes_t(_pattern_tip, _site_repeats, _rate_scalers, _cpu);
}

unsigned int test_case_t::attributes() const {
  return cpu_attributes() | misc_attributes();
}

unsigned int test_case_t::misc_attributes() const {
  return _pattern_tip ? PLL_ATTRIB_PATTERN_TIP
                      : 0 | _site_repeats
                            ? PLL_ATTRIB_SITE_REPEATS
                            : 0 | _rate_scalers ? PLL_ATTRIB_RATE_SCALERS : 0;
}

unsigned int test_case_t::cpu_attributes() const {
  if (test_cpu_t::none == _cpu) {
    return PLL_ATTRIB_ARCH_CPU;
  }
  if (test_cpu_t::sse == _cpu) {
    return PLL_ATTRIB_ARCH_SSE;
  }
  if (test_cpu_t::avx == _cpu) {
    return PLL_ATTRIB_ARCH_AVX;
  }
  if (test_cpu_t::avx2 == _cpu) {
    return PLL_ATTRIB_ARCH_AVX2;
  }
  if (test_cpu_t::avx512 == _cpu) {
    return PLL_ATTRIB_ARCH_AVX512;
  }
  throw std::runtime_error("Unrecognized CPU type");
}

} // namespace dks

std::ostream &operator<<(std::ostream &stream, const dks::test_cpu_t &cpu) {
  if (cpu == dks::test_cpu_t::none) {
    stream << "none";
  } else if (cpu == dks::test_cpu_t::sse) {
    stream << "sse";
  } else if (cpu == dks::test_cpu_t::avx) {
    stream << "avx";
  } else if (cpu == dks::test_cpu_t::avx2) {
    stream << "avx2";
  } else if (cpu == dks::test_cpu_t::avx512) {
    stream << "avx512";
  }
  return stream;
}

std::ostream &operator<<(std::ostream &stream,
                         const dks::attributes_t &attribs) {
  stream << "{\"cpu\": \"" << attribs.simd << "\"" << std::boolalpha
         << ", \"pattern tip\": \"" << attribs.pattern_tip << "\""
         << ", \"rate scalers\": \"" << attribs.rate_scalers << "\""
         << ", \"site repeats\": \"" << attribs.site_repeats << "\"}";
  return stream;
}
