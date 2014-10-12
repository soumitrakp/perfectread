
#ifndef __SEQUENCE_HANDLER_HPP__
#define __SEQUENCE_HANDLER_HPP__

#include "perfectread.hpp"

struct rule_zero : public filter_t {
  rule_zero(const params& par, const binary_query& bq, filter_t* prev = 0) : filter_t(par, bq, prev) { }
  bool satisfy(const qual_mer_dna& m) {
    //std::cout << m.to_str() << "=" << _query.check(m) << " vs " << _par.ce << "\n";
    size_t my_count=_query.check(m);
    return my_count >= _par.ce;
  }
};

struct rule_one : public filter_t {
  rule_one(const params& par, const binary_query& bq, filter_t* prev = 0) : filter_t(par, bq, prev) { }
  bool satisfy(const qual_mer_dna& m) {
    size_t my_count=_query.check(m);
    if (my_count < _par.cg) return false;
    for (const char *q=m.qual_str().c_str(); *q; ++q) {
      if (*q < _par.qg) return false;
    }
    return true;
  }
};



template<typename MerType>
class sequence_handler {
  sequence_parser&   _parser;
  const binary_query _query;
  const params&      _par;
  filter_t*          filter_;
  bool               canonical_;
  MerType            _mer; // mer
  MerType            _rev_mer; // reverse complement mer
  jflib::omstream & _output;
  jflib::omstream & _details;

public:
  typedef MerType mer_type;

  sequence_handler(sequence_parser& parser, const mapped_file& bm, const file_header& header, const params& par, jflib::omstream & output, jflib::omstream & details, bool canonical = false) :
    _parser(parser), _query(bm.base() + header.offset(), header.key_len(), header.counter_len(), header.matrix(),
        header.size() - 1, bm.length() - header.offset()), _par(par), canonical_(canonical), _output(output), _details(details) {

      rule_zero *rule0 = new rule_zero(par, _query);
      filter_ = new rule_one(par, _query, rule0);
    }

  int has_error(const char* seq, const char *qual, int n) {
    //std::cout << "seq::" << seq << "\n";
    //std::cout << "qul::" << qual << "\n";
    int k = qual_mer_dna::k();
    int filled_ = 0;
    int i=0;
    while (seq[i]) {
      int code = _mer.code(seq[i]);
      if(code < 0) {
        return 1;
      }
      char q = qual[i];
      _mer.shift_left(code, q);
      _rev_mer.shift_right(_rev_mer.complement(code), q);
      //std::cout << "mer:" << _mer.to_str() << ", qual:" << _mer.qual_str() << "\n";
      //std::cout << "rev:" << _rev_mer.to_str() << ", qual:" << _rev_mer.qual_str() << "\n";
      i++; filled_++;
      if (filled_>=k || seq[i] == 0) {
        //std::cout << "using mer:" << _mer.to_str() << "\n";
        qual_mer_dna &mer = (_mer < _rev_mer) ? _mer : _rev_mer;
        if(!(*filter_)(mer)) {
          return 1;
        }
        filled_-= k/2; //1; //k/8;
      }
    }
    return 0;
  }

  bool process(int thid) {
    //std::cout << "thread::" << thid << "\n";
    sequence_parser::job job(_parser);
    if (job.is_empty()) return false;
    for (size_t i=0; i<job->nb_filled; ++i) {
        const std::string& header = job->data[i].header;
        const std::string& seq    = job->data[i].seq;
        const std::string& qual   = job->data[i].qual;
        
      if (seq.size() != qual.size()) {
        std::cerr << "Error: Insufficient quality information. Quiting...\n";
        exit(0);
      }
      if (!has_error(seq.c_str(), qual.c_str(), seq.size())) {
       // std::cout << "READ IS PERFECT" << "\n";
        _output << "@" << header << "\n" << seq << "\n";
        _output << "+" << header << "\n" << qual << "\n";
        _output << jflib::endr;
      } else {
        //std::cout << "READ IS ERRONEOUS" << "\n";
      }
    }
    return true;
  }
};


#endif /* __SEQUENCE_HANDLER_HPP__ */
