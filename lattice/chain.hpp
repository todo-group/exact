// Copyright (C) 1997-2017 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef LATTICE_CHAIN_HPP
#define LATTICE_CHAIN_HPP

namespace lattice {

class chain {
public:
  chain(unsigned int L) : length_(L) {}
  unsigned int get_length() const { return length_; }
  unsigned int num_sites() const { return get_length(); }
  unsigned int num_bonds() const { return get_length(); }
  unsigned int source(unsigned int b) const { return b; }
  unsigned int target(unsigned int b) const { return (b + 1) % get_length(); }
  unsigned int num_neighbors() const { return 2; }
  unsigned int neighbor(unsigned int s, unsigned int k) const {
    return (get_length() + s + 1 - 2 * k) % get_length();
  }
  double site_phase(unsigned int s) const { return 2.0 * (s & 1) - 1.0; }
  double bond_phase(unsigned int b) const { return site_phase(b); }
private:
  unsigned int length_;
};

} // end namespace lattice

#endif // LATTICE_CHAIN_HPP
