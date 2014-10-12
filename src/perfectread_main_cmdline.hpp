/***** This code was generated by Yaggo. Do not edit ******/

/*  This file is part of Jellyfish.

    Jellyfish is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Jellyfish is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Jellyfish.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __ED_MAIN_CMDLINE_HPP__
#define __ED_MAIN_CMDLINE_HPP__

#include <stdint.h>
#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>
#include <errno.h>
#include <string.h>
#include <stdexcept>
#include <string>
#include <limits>
#include <vector>
#include <iostream>
#include <sstream>

class ed_main_cmdline {
 // Boiler plate stuff. Conversion from string to other formats
  static bool adjust_double_si_suffix(double &res, const char *suffix) {
    if(*suffix == '\0')
      return true;
    if(*(suffix + 1) != '\0')
      return false;

    switch(*suffix) {
    case 'a': res *= 1e-18; break;
    case 'f': res *= 1e-15; break;
    case 'p': res *= 1e-12; break;
    case 'n': res *= 1e-9;  break;
    case 'u': res *= 1e-6;  break;
    case 'm': res *= 1e-3;  break;
    case 'k': res *= 1e3;   break;
    case 'M': res *= 1e6;   break;
    case 'G': res *= 1e9;   break;
    case 'T': res *= 1e12;  break;
    case 'P': res *= 1e15;  break;
    case 'E': res *= 1e18;  break;
    default: return false;
    }
    return true;
  }

  static double conv_double(const char *str, ::std::string &err, bool si_suffix) {
    char *endptr = 0;
    errno = 0;
    double res = strtod(str, &endptr);
    if(errno) {
      err.assign(strerror(errno));
      return (double)0.0;
    }
    bool invalid =
      si_suffix ? !adjust_double_si_suffix(res, endptr) : *endptr != '\0';
    if(invalid) {
      err.assign("Invalid character");
      return (double)0.0;
    }
    return res;
  }

  static int conv_enum(const char* str, ::std::string& err, const char* const strs[]) {
    int res = 0;
    for(const char* const* cstr = strs; *cstr; ++cstr, ++res)
      if(!strcmp(*cstr, str))
        return res;
    err += "Invalid constant '";
    err += str;
    err += "'. Expected one of { ";
    for(const char* const* cstr = strs; *cstr; ++cstr) {
      if(cstr != strs)
        err += ", ";
      err += *cstr;
    }
    err += " }";
    return -1;
  }

  template<typename T>
  static bool adjust_int_si_suffix(T &res, const char *suffix) {
    if(*suffix == '\0')
      return true;
    if(*(suffix + 1) != '\0')
      return false;

    switch(*suffix) {
    case 'k': res *= (T)1000; break;
    case 'M': res *= (T)1000000; break;
    case 'G': res *= (T)1000000000; break;
    case 'T': res *= (T)1000000000000; break;
    case 'P': res *= (T)1000000000000000; break;
    case 'E': res *= (T)1000000000000000000; break;
    default: return false;
    }
    return true;
  }

  template<typename T>
  static T conv_int(const char *str, ::std::string &err, bool si_suffix) {
    char *endptr = 0;
    errno = 0;
    long long int res = strtoll(str, &endptr, 0);
    if(errno) {
      err.assign(strerror(errno));
      return (T)0;
    }
    bool invalid =
      si_suffix ? !adjust_int_si_suffix(res, endptr) : *endptr != '\0';
    if(invalid) {
      err.assign("Invalid character");
      return (T)0;
    }
    if(res > ::std::numeric_limits<T>::max() ||
       res < ::std::numeric_limits<T>::min()) {
      err.assign("Value out of range");
      return (T)0;
    }
    return (T)res;
  }

  template<typename T>
  static T conv_uint(const char *str, ::std::string &err, bool si_suffix) {
    char *endptr = 0;
    errno = 0;
    while(isspace(*str)) { ++str; }
    if(*str == '-') {
      err.assign("Negative value");
      return (T)0;
    }
    unsigned long long int res = strtoull(str, &endptr, 0);
    if(errno) {
      err.assign(strerror(errno));
      return (T)0;
    }
    bool invalid =
      si_suffix ? !adjust_int_si_suffix(res, endptr) : *endptr != '\0';
    if(invalid) {
      err.assign("Invalid character");
      return (T)0;
    }
    if(res > ::std::numeric_limits<T>::max() ||
       res < ::std::numeric_limits<T>::min()) {
      err.assign("Value out of range");
      return (T)0;
    }
    return (T)res;
  }

  template<typename T>
  static ::std::string vec_str(const std::vector<T> &vec) {
    ::std::ostringstream os;
    for(typename ::std::vector<T>::const_iterator it = vec.begin();
        it != vec.end(); ++it) {
      if(it != vec.begin())
        os << ",";
      os << *it;
    }
    return os.str();
  }

  class string : public ::std::string {
  public:
    string() : ::std::string() {}
    explicit string(const ::std::string &s) : std::string(s) {}
    explicit string(const char *s) : ::std::string(s) {}
    int as_enum(const char* const strs[]) {
      ::std::string err;
      int res = conv_enum((const char*)this->c_str(), err, strs);
      if(!err.empty())
        throw ::std::runtime_error(err);
      return res;
    }


    uint32_t as_uint32_suffix() const { return as_uint32(true); }
    uint32_t as_uint32(bool si_suffix = false) const {
      ::std::string err;
      uint32_t res = conv_uint<uint32_t>((const char*)this->c_str(), err, si_suffix);
      if(!err.empty()) {
        ::std::string msg("Invalid conversion of '");
        msg += *this;
        msg += "' to uint32_t: ";
        msg += err;
        throw ::std::runtime_error(msg);
      }
      return res;
    }
    uint64_t as_uint64_suffix() const { return as_uint64(true); }
    uint64_t as_uint64(bool si_suffix = false) const {
      ::std::string err;
      uint64_t res = conv_uint<uint64_t>((const char*)this->c_str(), err, si_suffix);
      if(!err.empty()) {
        ::std::string msg("Invalid conversion of '");
        msg += *this;
        msg += "' to uint64_t: ";
        msg += err;
        throw ::std::runtime_error(msg);
      }
      return res;
    }
    int32_t as_int32_suffix() const { return as_int32(true); }
    int32_t as_int32(bool si_suffix = false) const {
      ::std::string err;
      int32_t res = conv_int<int32_t>((const char*)this->c_str(), err, si_suffix);
      if(!err.empty()) {
        ::std::string msg("Invalid conversion of '");
        msg += *this;
        msg += "' to int32_t: ";
        msg += err;
        throw ::std::runtime_error(msg);
      }
      return res;
    }
    int64_t as_int64_suffix() const { return as_int64(true); }
    int64_t as_int64(bool si_suffix = false) const {
      ::std::string err;
      int64_t res = conv_int<int64_t>((const char*)this->c_str(), err, si_suffix);
      if(!err.empty()) {
        ::std::string msg("Invalid conversion of '");
        msg += *this;
        msg += "' to int64_t: ";
        msg += err;
        throw ::std::runtime_error(msg);
      }
      return res;
    }
    int as_int_suffix() const { return as_int(true); }
    int as_int(bool si_suffix = false) const {
      ::std::string err;
      int res = conv_int<int>((const char*)this->c_str(), err, si_suffix);
      if(!err.empty()) {
        ::std::string msg("Invalid conversion of '");
        msg += *this;
        msg += "' to int_t: ";
        msg += err;
        throw ::std::runtime_error(msg);
      }
      return res;
    }
    long as_long_suffix() const { return as_long(true); }
    long as_long(bool si_suffix = false) const {
      ::std::string err;
      long res = conv_int<long>((const char*)this->c_str(), err, si_suffix);
      if(!err.empty()) {
        ::std::string msg("Invalid conversion of '");
        msg += *this;
        msg += "' to long_t: ";
        msg += err;
        throw ::std::runtime_error(msg);
      }
      return res;
    }
    double as_double_suffix() const { return as_double(true); }
    double as_double(bool si_suffix = false) const {
      ::std::string err;
      double res = conv_double((const char*)this->c_str(), err, si_suffix);
      if(!err.empty()) {
        ::std::string msg("Invalid conversion of '");
        msg += *this;
        msg += "' to double_t: ";
        msg += err;
        throw ::std::runtime_error(msg);
      }
      return res;
    }
  };

public:
  uint64_t                       size_arg;
  bool                           size_given;
  uint32_t                       threads_arg;
  bool                           threads_given;
  uint32_t                       Files_arg;
  bool                           Files_given;
  const char *                   generator_arg;
  bool                           generator_given;
  uint32_t                       Generators_arg;
  bool                           Generators_given;
  const char *                   shell_arg;
  bool                           shell_given;
  const char *                   output_arg;
  bool                           output_given;
  const char *                   db_arg;
  bool                           db_given;
  uint32_t                       good_qual_arg;
  bool                           good_qual_given;
  uint32_t                       excel_qual_arg;
  bool                           excel_qual_given;
  uint32_t                       good_count_arg;
  bool                           good_count_given;
  uint32_t                       excel_count_arg;
  bool                           excel_count_given;
  const char *                   timing_arg;
  bool                           timing_given;
  ::std::vector<const char *>    file_arg;
  typedef ::std::vector<const char *>::iterator file_arg_it;
  typedef ::std::vector<const char *>::const_iterator file_arg_const_it;

  enum {
    START_OPT = 1000,
    FULL_HELP_OPT,
    USAGE_OPT,
    TIMING_OPT
  };

  ed_main_cmdline() :
    size_arg(), size_given(false),
    threads_arg(1), threads_given(false),
    Files_arg(1), Files_given(false),
    generator_arg(""), generator_given(false),
    Generators_arg(1), Generators_given(false),
    shell_arg(""), shell_given(false),
    output_arg("mer_counts.jf"), output_given(false),
    db_arg("mer_counts.jf"), db_given(false),
    good_qual_arg(45), good_qual_given(false),
    excel_qual_arg(73), excel_qual_given(false),
    good_count_arg(1), good_count_given(false),
    excel_count_arg(1), excel_count_given(false),
    timing_arg(""), timing_given(false),
    file_arg()
  { }

  ed_main_cmdline(int argc, char* argv[]) :
    size_arg(), size_given(false),
    threads_arg(1), threads_given(false),
    Files_arg(1), Files_given(false),
    generator_arg(""), generator_given(false),
    Generators_arg(1), Generators_given(false),
    shell_arg(""), shell_given(false),
    output_arg("mer_counts.jf"), output_given(false),
    db_arg("mer_counts.jf"), db_given(false),
    good_qual_arg(45), good_qual_given(false),
    excel_qual_arg(73), excel_qual_given(false),
    good_count_arg(1), good_count_given(false),
    excel_count_arg(1), excel_count_given(false),
    timing_arg(""), timing_given(false),
    file_arg()
  { parse(argc, argv); }

  void parse(int argc, char* argv[]) {
    static struct option long_options[] = {
      {"size", 1, 0, 's'},
      {"threads", 1, 0, 't'},
      {"Files", 1, 0, 'F'},
      {"generator", 1, 0, 'g'},
      {"Generators", 1, 0, 'G'},
      {"shell", 1, 0, 'S'},
      {"output", 1, 0, 'o'},
      {"db", 1, 0, 'd'},
      {"good-qual", 1, 0, 'Q'},
      {"excellent-qual", 1, 0, 'q'},
      {"good-count", 1, 0, 'C'},
      {"excellent-count", 1, 0, 'c'},
      {"timing", 1, 0, TIMING_OPT},
      {"help", 0, 0, 'h'},
      {"full-help", 0, 0, FULL_HELP_OPT},
      {"usage", 0, 0, USAGE_OPT},
      {"version", 0, 0, 'V'},
      {0, 0, 0, 0}
    };
    static const char *short_options = "hVs:t:Q:q:C:c:F:g:G:S:d:o:";

    ::std::string err;
#define CHECK_ERR(type,val,which) if(!err.empty()) { ::std::cerr << "Invalid " #type " '" << val << "' for [" which "]: " << err << "\n"; exit(1); }
    while(true) {
      int index = -1;
      int c = getopt_long(argc, argv, short_options, long_options, &index);
      if(c == -1) break;
      switch(c) {
      case ':':
        ::std::cerr << "Missing required argument for "
                  << (index == -1 ? ::std::string(1, (char)optopt) : std::string(long_options[index].name))
                  << ::std::endl;
        exit(1);
      case 'h':
        ::std::cout << usage() << "\n\n" << help() << std::endl;
        exit(0);
      case USAGE_OPT:
        ::std::cout << usage() << "\nUse --help for more information." << std::endl;
        exit(0);
      case 'V':
        print_version();
        exit(0);
      case '?':
        ::std::cerr << "Use --usage or --help for some help\n";
        exit(1);
      case FULL_HELP_OPT:
        ::std::cout << usage() << "\n\n" << help() << "\n\n" << hidden() << std::endl;
        exit(0);
      case 's':
        size_given = true;
        size_arg = conv_uint<uint64_t>((const char*)optarg, err, true);
        CHECK_ERR(uint64_t, optarg, "-s, --size=uint64")
        break;
      case 't':
        threads_given = true;
        threads_arg = conv_uint<uint32_t>((const char*)optarg, err, false);
        CHECK_ERR(uint32_t, optarg, "-t, --threads=uint32")
        break;
      case 'F':
        Files_given = true;
        Files_arg = conv_uint<uint32_t>((const char*)optarg, err, false);
        CHECK_ERR(uint32_t, optarg, "-F, --Files=uint32")
        break;
      case 'g':
        generator_given = true;
        generator_arg = optarg;
        break;
      case 'G':
        Generators_given = true;
        Generators_arg = conv_uint<uint32_t>((const char*)optarg, err, false);
        CHECK_ERR(uint32_t, optarg, "-G, --Generators=uint32")
        break;
      case 'S':
        shell_given = true;
        shell_arg = optarg;
     case 'o':
        output_given = true;
        output_arg = optarg;
        break;
      case 'd':
        db_given = true;
        db_arg = optarg;
        break;
      case 'Q':
        good_qual_given = true;
        good_qual_arg = conv_uint<uint32_t>((const char*)optarg, err, false);
        CHECK_ERR(uint8_t, optarg, "-Q, --good-qual=uint8")
        break;
      case 'q':
        excel_qual_given = true;
        excel_qual_arg = conv_uint<uint32_t>((const char*)optarg, err, false);
        CHECK_ERR(uint8_t, optarg, "-q, --excellent-qual=uint8")
        break;
      case 'C':
        good_count_given = true;
        good_count_arg = conv_uint<uint32_t>((const char*)optarg, err, false);
        CHECK_ERR(uint32_t, optarg, "-C, --good-count=uint32")
        break;
      case 'c':
        excel_count_given = true;
        excel_count_arg = conv_uint<uint32_t>((const char*)optarg, err, false);
        CHECK_ERR(uint32_t, optarg, "-c, --excellent-count=uint32")
        break;
      case TIMING_OPT:
        timing_given = true;
        timing_arg = optarg;
        break;
      }
    }

    // Check that required switches are present
    if(!size_given)
      error("[-s, --size=uint64] required switch");

    if(!db_given)
      error("[-d, --db=string] required switch");

    // Parse arguments
    if(argc - optind < 0)
      error("Requires at least 0 argument.");
    for( ; optind < argc; ++optind) {
      file_arg.push_back(argv[optind]);
    }
  }

#define ed_main_cmdline_USAGE "Usage: perfectread [options] file:path+"

  const char * usage() const { return ed_main_cmdline_USAGE; }
  void error(const char *msg) {
    ::std::cerr << "Error: " << msg << "\n" << usage()
              << "\nUse --help for more information"
              << ::std::endl;
    exit(1);
  }

#define ed_main_cmdline_HELP "detect erroneous reads in fasta or fastq files\n\n" \
  "Options (default value in (), *required):\n" \
  " -s, --size=uint64                       *Number of reads in fasta\n" \
  " -d, --db=string                         *Input kmer database file created using count subcommand\n" \
  " -F, --Files=uint32                       Number files open simultaneously (1)\n" \
  " -g, --generator=path                     File of commands generating fast[aq]\n" \
  " -G, --Generators=uint32                  Number of generators run simultaneously (1)\n" \
  " -S, --shell=string                       Shell used to run generator commands ($SHELL or /bin/sh)\n" \
  " -t, --threads=uint32                     Number of threads (1)\n" \
  " -o, --output=string                      Output file (mer_counts.jf)\n" \
  " -Q, --good-qual=uint32                   Good quality value in integer (45)\n" \
  " -q, --excellent-qual=uint32              Excellent quality value in integer (73)\n" \
  " -C, --good-count=uint32                  Good count in integer (1)\n" \
  " -c, --excellent-count=uint32             Excellent count in integer (1)\n" \
  "     --timing=Timing file                 Print timing information\n" \
  "     --usage                              Usage\n" \
  " -h, --help                               This message\n" \
  "     --full-help                          Detailed help\n" \
  " -V, --version                            Version"
  const char * help() const { return ed_main_cmdline_HELP; }

#define ed_main_cmdline_HIDDEN "Hidden options:\n" \
  "     --no-merge                           Do not merge files intermediary files (false)\n" \
  "     --no-unlink                          Do not unlink intermediary files after automatic merging (false)\n" \
  "     --no-write                           Don't write database (false)"

  const char * hidden() const { return ed_main_cmdline_HIDDEN; }
  void print_version(::std::ostream &os = std::cout) const {
#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.0.0"
#endif
    os << PACKAGE_VERSION << "\n";
  }
  void dump(::std::ostream &os = std::cout) {
    os << "size_given:" << size_given << " size_arg:" << size_arg << "\n";
    os << "threads_given:" << threads_given << " threads_arg:" << threads_arg << "\n";
    os << "Files_given:" << Files_given << " Files_arg:" << Files_arg << "\n";
    os << "generator_given:" << generator_given << " generator_arg:" << generator_arg << "\n";
    os << "Generators_given:" << Generators_given << " Generators_arg:" << Generators_arg << "\n";
    os << "shell_given:" << shell_given << " shell_arg:" << shell_arg << "\n";
    os << "output_given:" << output_given << " output_arg:" << output_arg << "\n";
    os << "db_given:" << db_given << " db_arg:" << db_arg << "\n";
    os << "good_qual_given:" << good_qual_given << " good_qual_arg:" << good_qual_arg << "\n";
    os << "excel_qual_given:" << excel_qual_given << " excel_qual_arg:" << excel_qual_arg << "\n";
    os << "good_count_given:" << good_count_given << " good_count_arg:" << good_count_arg << "\n";
    os << "excel_count_given:" << excel_count_given << " excel_count_arg:" << excel_count_arg << "\n";
    os << "timing_given:" << timing_given << " timing_arg:" << timing_arg << "\n";
    os << "file_arg:" << vec_str(file_arg) << "\n";
  }
};
#endif // __ED_MAIN_CMDLINE_HPP__"
