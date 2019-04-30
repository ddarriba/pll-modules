#pragma once
#include "test_case.h"
#include <unordered_map>

namespace dks {
typedef std::unordered_map<test_kernel_t, double> kernel_weight_t;
typedef std::unordered_map<attributes_t, benchmark_time_t> attributes_time_t;

kernel_weight_t suggest_weights(const msa_t& msa);
kernel_weight_t suggest_weights_2(const msa_t& msa);

attributes_time_t select_kernel_fast_verbose(const model_t &model, const msa_t &msa,
                           const kernel_weight_t &);
attributes_time_t select_kernel_slow_verbose(const model_t &model, const msa_t &msa,
                           const kernel_weight_t &);
attributes_t select_kernel(const pll_partition_t *, const pll_msa_t *,
                           const kernel_weight_t &, bool fast);

attributes_t select_kernel(const model_t &, const msa_t &,
                           const kernel_weight_t &, bool fast);

attributes_time_t select_kernel_verbose(const model_t &, const msa_t &,
                                        const kernel_weight_t &, bool fast);
} // namespace dks
