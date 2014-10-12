#include <chrono>

using std::chrono::system_clock;
using std::chrono::duration;
using std::chrono::duration_cast;
template<typename DtnType>
inline double as_seconds(DtnType dtn) { return duration_cast<duration<double>>(dtn).count(); }

#include "perfectread.hpp"
#include "sequence_handler.hpp"
#include "perfectread_main_cmdline.hpp"

typedef sequence_handler<qual_mer_dna> seq_handler;

static ed_main_cmdline args;

template<typename PathIterator>
class error_detect : public jellyfish::thread_exec {
  jellyfish::stream_manager<PathIterator> _streams;
  sequence_parser       _parser;
  mapped_file&          _map;
  file_header&		      _header;
  params& 			        _par;
  std::string           _prefix;
  jflib::o_multiplexer* _output;
  jflib::o_multiplexer* _log;

private:
  // Open the data (error corrected reads) and log files. Default to
  // STDOUT and STDERR if none specified.
  std::ostream* open_file(const std::string prefix, const char* suffix,
      const std::string def) {
    std::ostream* res;
    std::string file;
    if(prefix.empty())
      file = def;
    else {
      file = prefix;
      file += suffix;
    }
    res = new std::ofstream(file.c_str());
    if(!res->good()) {
      std::cerr << "Failed to open file '" << file << "'\n";
      exit(-1);
    }
    res->exceptions(std::ios::eofbit|std::ios::failbit|std::ios::badbit);
    return res;
  }

public:
  error_detect(int nb_threads, uint32_t nb_sequences, mapped_file& bm, file_header& header,
      PathIterator file_begin, PathIterator file_end, params& par) :
    _streams(file_begin, file_end),
    _parser(3 * nb_threads, nb_sequences, _streams.nb_streams(),  _streams),
    _map(bm), _header(header), _par(par) {
      std::string first_file(*file_begin);
      _prefix = first_file.substr(0, first_file.find_last_of(".")) + "-perfect";
    }

  void do_it(int nb_threads) {
    std::unique_ptr<std::ostream> out(open_file(_prefix, ".fastq", "/dev/fd/1"));
    std::unique_ptr<std::ostream> log(open_file(_prefix, ".log", "/dev/fd/2"));
    std::unique_ptr<jflib::o_multiplexer> out_m(new jflib::o_multiplexer(out.get(), 3 * nb_threads, 1024));
    std::unique_ptr<jflib::o_multiplexer> log_m(new jflib::o_multiplexer(log.get(), 3 * nb_threads, 1024));
    _output = out_m.get();
    _log = log_m.get();
    exec_join(nb_threads);
  }

  virtual void start(int thid) {
    jflib::omstream op(output());
    jflib::omstream lg(log());
    seq_handler handler(_parser, _map, _header, _par, op, lg);
    while (handler.process(thid));
    op.close();
    lg.close();
  }

  jflib::o_multiplexer& output() { return *_output; }
  jflib::o_multiplexer& log() { return *_log; }

};

int main(int argc, char *argv[])
{
  auto start_time = system_clock::now();

  args.parse(argc, argv);

  ofstream_default out(args.output_given ? args.output_arg : 0, std::cout);
  if(!out.good())
    std::cerr << "Error opening output file '" << args.output_arg << "'";

  std::ifstream in(args.db_arg, std::ios::in|std::ios::binary);
  jellyfish::file_header header(in);
  if(!in.good())
    std::cerr << "Failed to parse header of file '" << args.db_arg << "'";
  qual_mer_dna::k(header.key_len() / 2);
  std::cout << "Reading from database ... (k=" << qual_mer_dna::k() << ")\n";
  if(header.format() != binary_dumper::format) {
    std::cerr << "Unsupported format '" << header.format() << "'. Must be a bloom counter or binary list.";
  }
  jellyfish::mapped_file binary_map(args.db_arg);

  std::cout << "loaded database \n";

  auto after_init_time = system_clock::now();

  params par(args.good_qual_arg, args.excel_qual_arg, args.good_count_arg, args.excel_count_arg);

  error_detect<file_vector::const_iterator>  ed(args.threads_arg, args.size_arg, binary_map, header, 
      args.file_arg.begin(), args.file_arg.end(), par);

  ed.do_it(args.threads_arg);

  auto after_count_time = system_clock::now();

  std::cout << "Init     " << as_seconds(after_init_time - start_time) << "\n"
      << "Correction " << as_seconds(after_count_time - after_init_time) << "\n";

  if(args.timing_given) {
    std::ofstream timing_file(args.timing_arg);
    timing_file << "Init     " << as_seconds(after_init_time - start_time) << "\n"
      << "Correction " << as_seconds(after_count_time - after_init_time) << "\n";
  }

  return 0;
}
