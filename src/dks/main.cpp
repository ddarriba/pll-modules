#include "benchmark.h"
#include <algorithm>
#include <getopt.h>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <unistd.h>

using namespace std;

option cli_options[] = {
    {"msa", required_argument, 0, 0},
    {"states", required_argument, 0, 0},
    {0, 0, 0, 0},
};

int main(int argc, char **argv) {
  int opt_index;
  string filename;
  char c;
  size_t states = 4;
  while ((c = getopt_long(argc, argv, "", cli_options, &opt_index)) == 0) {
    switch (opt_index) {
    case 0:
      filename = string(optarg);
      break;
    case 1:
      states = std::stoi(optarg);
      break;
    }
  }
  if (filename.empty()) {
    cout << "Please give an msa file with the --msa switch" << endl;
    return 1;
  }

  dks::msa_t msa(filename, states);
  dks::model_t model(msa);

  dks::kernel_weight_t kw = dks::suggest_weights_2(msa);

  auto results = dks::select_kernel_verbose(model, msa, kw, false);

  vector<pair<dks::attributes_t, dks::benchmark_time_t>> sorted_times{
      results.begin(), results.end()};
  std::sort(sorted_times.begin(), sorted_times.end(),
            [](const decltype(sorted_times)::value_type &a,
               const decltype(sorted_times)::value_type &b) {
              return a.second < b.second;
            });

  for (const auto &kv : sorted_times) {
    cout << "{\"attributes\":" << kv.first
         << ", \"score\": " << kv.second.count() << "}" << endl;
  }

  return 0;
}
