// Copyright (C) 1997-2017 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef LATTICE_SQUARE_HPP
#define LATTICE_SQUARE_HPP

#include <boost/array.hpp>
#include <boost/tuple/tuple.hpp>
#include <vector>

namespace lattice {

class square {
public:
  square(unsigned int L) : length_x_(L), length_y_(L) { init(); }
  square(unsigned int Lx, unsigned int Ly) : length_x_(Lx), length_y_(Ly) { init(); }
  void init() {
    neighbors_.resize(num_sites());
    source_.resize(num_bonds());
    target_.resize(num_bonds());
    site_phase_.resize(num_sites());
    bond_phase_.resize(num_bonds());
    for (unsigned int s = 0; s < num_sites(); ++s) {
      int x, y;
      boost::tie(x, y) = index2xy(s);
      site_phase_[s] = 2.0 * ((x + y) % 2) - 1.0;
      for (unsigned int k = 0; k < 4; ++k) {
        int d = 1- int(k & 2);
        neighbors_[s][k] = ((k & 1) == 0) ? xy2index(x + d, y) : xy2index(x, y + d);
      }
    }
    for (unsigned int b = 0; b < num_bonds(); ++b) {
      unsigned int s = b / 2;
      int x, y;
      boost::tie(x, y) = index2xy(s);
      unsigned int t;
      if (b % 2 == 0) {
        t = xy2index(x + 1, y); // target right
        bond_phase_[b] = 2.0 * ((b / 2) % 2) - 1.0;
      } else {
        t = xy2index(x, y + 1); // target below
        bond_phase_[b] = 2.0 * ((b / length_x_ / 2) % 2) - 1.0;
      }
      source_[b] = s;
      target_[b] = t;
    }
  }
  unsigned int get_length_x() const { return length_x_; }
  unsigned int get_length_y() const { return length_y_; }
  unsigned int num_sites() const { return length_x_ * length_y_; }
  unsigned int num_bonds() const { return 2 * num_sites(); }
  unsigned int source(unsigned int b) const { return source_[b]; }
  unsigned int target(unsigned int b) const { return target_[b]; }
  unsigned int num_neighbors() const { return 4; }
  unsigned int neighbor(unsigned int s, unsigned int k) const { return neighbors_[s][k]; }
  double site_phase(unsigned int s) const { return site_phase_[s]; }
  double bond_phase(unsigned int b) const { return bond_phase_[b]; }
protected:
  boost::tuple<int, int> index2xy(unsigned int s) const {
    return boost::make_tuple(s % length_x_, s / length_x_);
  }
  unsigned int xy2index(int x, int y) const {
    x += length_x_;
    y += length_y_;
    return x % length_x_ + (y % length_y_) * length_x_;
  }
private:
  unsigned int length_x_, length_y_;
  std::vector<boost::array<unsigned int, 4> > neighbors_;
  std::vector<unsigned int> source_;
  std::vector<unsigned int> target_;
  std::vector<double> site_phase_;
  std::vector<double> bond_phase_;
};

} // end namespace lattice

#endif // LATTICE_SQUARE_HPP
