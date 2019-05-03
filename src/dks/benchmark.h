#pragma once
#include "pll_msa.h"
#include "test_case.h"
#include <unordered_map>

namespace dks {
typedef std::unordered_map<test_kernel_t, double> kernel_weight_t;
typedef std::unordered_map<attributes_t, benchmark_time_t> attributes_time_t;

typedef std::vector<std::vector<char>> msa_t;

kernel_weight_t suggest_weights(const msa_t &msa);
kernel_weight_t suggest_weights(const pllmod_msa_stats_t *msa, int sites,
                                int taxa);

template <typename T>
attributes_time_t select_kernel_verbose(const model_t &model, const T &msa,
                                        const msa_weight_t &weights,
                                        const pll_state_t *charmap,
                                        const kernel_weight_t &kw,
                                        attributes_generator_t att_gen);

attributes_t select_kernel(const pll_partition_t *, const pll_msa_t *,
                           const kernel_weight_t &);

attributes_t select_kernel(const model_t &, const msa_t &,
                           const kernel_weight_t &);

attributes_t select_kernel_auto(const pll_partition_t *pll_partition,
                                const pll_msa_t *pll_msa);

unsigned int select_kernel_auto(const msa_t &msa, const msa_weight_t &weights,
                                const pll_state_t *charmap, unsigned int states,
                                unsigned int rate_cats,
                                attributes_generator_t gen);

unsigned int select_kernel_auto(const msa_t &msa, const msa_weight_t &weights,
                                const pll_state_t *charmap, unsigned int states,
                                unsigned int rate_cats);

unsigned int select_kernel_auto(const std::vector<std::vector<double>> &clvs,
                                const msa_weight_t &weights,
                                unsigned int states, unsigned int rate_cats,
                                attributes_generator_t gen);

unsigned int select_kernel_auto(const std::vector<std::vector<double>> &clvs,
                                const msa_weight_t &weights,
                                unsigned int states, unsigned int rate_cats);

} // namespace dks
