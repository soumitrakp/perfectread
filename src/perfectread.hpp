#ifndef __PERFECT_READ_HPP__
#define __PERFECT_READ_HPP__

#include <cstdlib>
#include <unistd.h>
#include <assert.h>
#include <signal.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <memory>

#include <jellyfish/err.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/locks_pthread.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/whole_sequence_parser.hpp>
#include <jellyfish/generator_manager.hpp>
#include <jellyfish/fstream_default.hpp>
#include <jellyfish/hash_counter.hpp>
#include <jellyfish/binary_dumper.hpp>
#include <jellyfish/mapped_file.hpp>
#include <jflib/multiplexed_io.hpp>

#include "qual_mer_dna.hpp"
typedef jellyfish::cooperative::hash_counter<qual_mer_dna> mer_hash;
typedef jellyfish::binary_reader<qual_mer_dna, uint64_t> binary_reader;
typedef jellyfish::binary_query_base<qual_mer_dna, uint64_t> binary_query;
typedef jellyfish::binary_dumper<mer_hash::array> binary_dumper;

typedef std::vector<const char*> file_vector;

typedef jellyfish::whole_sequence_parser<jellyfish::stream_manager<file_vector::const_iterator> > sequence_parser;
typedef jellyfish::mapped_file mapped_file;
typedef jellyfish::file_header file_header;

struct params {
  char qg;
  char qe;
  size_t cg;
  size_t ce;
  params(size_t _qg, size_t _qe, size_t _cg, size_t _ce) : 
    qg(_qg), qe(_qe), cg(_cg), ce(_ce) { }
  ~params() {}
};


// k-mer filters. Organized in a linked list, interpreted as a &&
// (logical and). I.e. all filter must return true for the result to
// be true. By default, filter returns true.
struct filter_t {
  const params&         _par;
  const binary_query&   _query;
  filter_t* prev_;
  filter_t(const params& par, const binary_query& bq, filter_t* prev = 0) : _par(par), _query(bq), prev_(prev) { }
  virtual ~filter_t() { }
  virtual bool operator()(const qual_mer_dna& x) {
    if (prev_ && (*prev_)(x)) return true;
    return satisfy(x);
  }
  virtual bool satisfy(const qual_mer_dna& x) {
    return false;
  }
};

#endif 

