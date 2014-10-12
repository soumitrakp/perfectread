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

#ifndef __QUAL_MER_DNA_HPP__
#define __QUAL_MER_DNA_HPP__

#include <jellyfish/mer_dna.hpp>

template<typename T = uint64_t>
class qual_mer_base : public jellyfish::mer_dna_ns::mer_base<qual_mer_base<T> > {
public:
  typedef T base_type;
  typedef jellyfish::mer_dna_ns::mer_base<qual_mer_base<T> > super;

  qual_mer_base() : super(k_), _qual(new char[k_+1]) {
    memset(_qual, '!', (k_) * sizeof(char));
    _qual[k_] = '\0';
  }

  explicit qual_mer_base(unsigned int k) : super(k), _qual(new char[k+1]) {
    if(k != k_)
      throw std::length_error(jellyfish::mer_dna_ns::error_different_k);
    memset(_qual, '!', (k) * sizeof(char));
    _qual[k] = '\0';
  }

  qual_mer_base(const qual_mer_base &m) : super(m), _qual(new char[static_cast<const T*>(&m)->k()+1])
  {
    memcpy(_qual, m._qual, (static_cast<const T*>(&m)->k()+1) * sizeof(char));
  }

  template<typename U>
  qual_mer_base(const unsigned int k, const U& rhs) : super(k, rhs), _qual(new char[k+1]) {
    memset(_qual, '!', (k_) * sizeof(char));
    _qual[k_] = '\0';
  }

 // qual_mer_base& operator=(const char* s) { return super::operator=(s); }
 // qual_mer_base& operator=(const std::string& s) { return super::operator=(s); }

  static unsigned int k(); // { return k_; }
  static unsigned int k(unsigned int k) { std::swap(k, k_); return k; }


  ~qual_mer_base() {
    delete [] _qual;
  }

  char qual(unsigned int i) const { return _qual[i]; }

  T& operator=(const qual_mer_base& rhs) {
    super::operator=(rhs);
    memcpy(_qual, rhs._qual, super::k() * sizeof(char));
    return *static_cast<T*>(this);
  }

  base_type shift_left(int c, char q) {
    unsigned int last = super::k()-1;
    for (unsigned int i=0; i<last; i++) {
      _qual[i] = _qual[i+1];
    }
    _qual[last] = q;
    return super::shift_left(c);
  }

  base_type shift_right(int c, char q) {
    for (unsigned int i=super::k()-1; i>0; --i) {
      _qual[i] = _qual[i-1];
    }
    _qual[0] = q;
    return super::shift_right(c);
  }

  char shift_left(char c, char q) {
    int x = super::code(c);
    if(x == -1)
      return 'N';
    return rev_code(shift_left(x, q));
  }

  char shift_right(char c, char q) {
    int x = super::code(c);
    if(x == -1)
      return 'N';
    return rev_code(shift_right(x, q));
  }

  void reverse_complement() {
    super::reverse_complement();
    unsigned int k = super::k();
    for (unsigned int i=0; i<k/2; i++) {
      char t = _qual[i];
      _qual[i] = _qual[k-1-i];
      _qual[k-1-i] = t;
    }
  }

  std::string qual_str() const {
    std::string res(_qual);
    return res;
  }

protected:
  char *              _qual;
  static unsigned int k_;
};

template<typename T>
unsigned int qual_mer_base<T>::k_ = 22;
template<typename T>
unsigned int qual_mer_base<T>::k() { return k_; }

namespace jellyfish { namespace mer_dna_ns {
template<typename T>
struct mer_dna_traits<qual_mer_base<T> > {
  typedef T base_type;
};
} }

typedef qual_mer_base<uint64_t> qual_mer_dna;

#endif
