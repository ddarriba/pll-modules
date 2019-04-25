#pragma once
#include <iterator>
#include <pll.h>
#include <string>
#include <vector>

namespace dks {
typedef std::vector<char> sequence_t;
typedef std::string label_t;

class msa_t {
public:
  msa_t() = default;
  msa_t(const pll_msa_t *);
  msa_t(const pll_msa_t *, size_t);
  msa_t(const std::string & f) : msa_t(f, 4) {};
  msa_t(const std::string &, size_t);
  void init(const pll_msa_t *msa);
  size_t count() const;
  size_t length() const;
  const char *label(size_t i) const;
  const char *sequence(size_t i) const;
  const pll_state_t *char_map() const;
  bool valid() const;
  size_t states() const;
  void set_states(size_t);

  double column_entropy() const;
  double row_entropy() const;

protected:
  std::vector<sequence_t> _sequences;
  std::vector<label_t> _labels;
  size_t _states;
};

class msa_compressed_t : public msa_t {
public:
  msa_compressed_t(const msa_t &);
  const unsigned int *weights() const;

private:
  std::vector<unsigned int> _weights;
};
} // namespace dks
