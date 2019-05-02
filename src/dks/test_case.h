#pragma once
#include "model.h"
#include "msa.h"
#include "partition.h"
#include "tree.h"
#include <chrono>
#include <memory>
#include <pll.h>
#include <unordered_map>

#define DKS_ATTRIB_MASK 0x3ff

namespace dks {
enum test_cpu_t {
  none,
  sse,
  avx,
  avx2,
  avx512,
  invalid,
};

test_cpu_t test_cpu_from_attribs(uint32_t attribs);

enum test_kernel_t {
  partial,
  likelihood,
  derivative,
  pmatrix,
};

typedef std::chrono::duration<double> benchmark_time_t;

typedef std::unordered_map<test_kernel_t, benchmark_time_t> benchmark_result_t;

struct attributes_t {
  bool pattern_tip;
  bool site_repeats;
  bool rate_scalers;
  test_cpu_t simd;

  attributes_t() = default;

  attributes_t(bool pt, bool sr, bool rs, test_cpu_t s)
      : pattern_tip{pt}, site_repeats{sr}, rate_scalers{rs}, simd{s} {};

  bool operator==(const attributes_t &other) const {
    return pattern_tip == other.pattern_tip &&
           site_repeats == other.site_repeats &&
           rate_scalers == other.rate_scalers && simd == other.simd;
  }

  bool operator!=(const attributes_t &other) const { return !(*this == other); }

  int cpu_attrib() const {
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
  }

  int siterepeat_attrib() const {
    return site_repeats ? PLL_ATTRIB_SITE_REPEATS : 0;
  }

  int patterntip_attrib() const {
    return pattern_tip ? PLL_ATTRIB_PATTERN_TIP : 0;
  }

  int pll_attributes() const {
    return cpu_attrib() | siterepeat_attrib() | patterntip_attrib();
  }
};

class attributes_generator_t {
public:
  attributes_generator_t()
      : _off_flags{PLL_ATTRIB_AB_MASK | PLL_ATTRIB_AB_FLAG |
                   PLL_ATTRIB_RATE_SCALERS | PLL_ATTRIB_ARCH_AVX512},
        _on_flags{0}, _max{(1 << 11) - 1}, _state{0} {};

  attributes_t next() {
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

  attributes_t end() {
    return attributes_t(false, false, false, dks::test_cpu_t::invalid);
  }

  void disable(int attribs) { _off_flags |= attribs; }
  void enable(int attribs) { _on_flags |= attribs; }

private:
  inline bool check_cases() const {
    // do some bit math to check that a bit is only set if the corrisponding
    // flag is set

    return ((_off_flags & _state) | (_on_flags & ~_state)) & DKS_ATTRIB_MASK;
  }

  bool check_xor(int attrib) const {
    return __builtin_popcount(attrib & _state) <= 1;
  }

  bool valid() const {
    if (!check_xor(PLL_ATTRIB_SITE_REPEATS | PLL_ATTRIB_PATTERN_TIP)) {
      return false;
    }
    if (!check_xor(PLL_ATTRIB_ARCH_MASK)) {
      return false;
    }
    return !check_cases();
  }

  uint32_t _off_flags;
  uint32_t _on_flags;
  uint32_t _max;
  uint32_t _state;
};

class test_case_t {
public:
  test_case_t()
      : _cpu{test_cpu_t::none}, _trials{30}, _random_seed{0},
        _pattern_tip{false}, _site_repeats{false}, _rate_scalers{false} {}

  test_case_t(test_cpu_t cpu, bool pt, bool sr, bool rs, uint64_t seed)
      : _cpu{cpu}, _trials{30}, _random_seed{seed}, _pattern_tip{pt},
        _site_repeats{sr}, _rate_scalers{rs} {}

  test_case_t(test_cpu_t cpu) : test_case_t{cpu, 0, 0, 0, 0} {}

  test_case_t(const attributes_t &attribs)
      : test_case_t(attribs.simd, attribs.pattern_tip, attribs.site_repeats,
                    attribs.rate_scalers, 0){};

  benchmark_result_t benchmark(const msa_t &, const model_t &);
  benchmark_time_t benchmark_partials(partition_t &partition,
                                      const model_t &model);
  benchmark_time_t benchmark_likelihood(partition_t &partition,
                                        const model_t &model);
  benchmark_time_t benchmark_pmatrix(partition_t &partition,
                                     const model_t &model);
  benchmark_time_t benchmark_derivative(partition_t &partition,
                                        const model_t &model);
  benchmark_time_t benchmark_update_site_repeats(partition_t &partition,
                                                 const model_t &model);
  attributes_t attributes_struct() const;
  unsigned int attributes() const;
  unsigned int cpu_attributes() const;
  unsigned int misc_attributes() const;

private:
  test_cpu_t _cpu;
  size_t _trials;
  uint64_t _random_seed;
  bool _pattern_tip;
  bool _site_repeats;
  bool _rate_scalers;
};
} // namespace dks

namespace std {
template <> struct hash<dks::attributes_t> {
  typedef dks::attributes_t argument_type;
  typedef size_t result_type;
  result_type operator()(const argument_type &s) const noexcept {
    return (s.pattern_tip << 0) ^ (s.site_repeats << 1) ^
           (s.rate_scalers << 2) ^ (s.simd << 3);
  }
};
template <> struct hash<dks::test_kernel_t> {
  typedef dks::test_kernel_t argument_type;
  typedef size_t result_type;
  result_type operator()(const argument_type &s) const noexcept {
    return static_cast<size_t>(s);
  }
};
} // namespace std

std::ostream &operator<<(std::ostream &stream, const dks::test_cpu_t &cpu);
std::ostream &operator<<(std::ostream &stream,
                         const dks::attributes_t &attribs);
