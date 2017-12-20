// Copyright (C) 1997-2017 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef LATTICE_CUBIC_HPP
#define LATTICE_CUBIC_HPP

#include <vector>
#include <boost/array.hpp>
#include <boost/tuple/tuple.hpp>

namespace lattice {

class cubic {
public:
  cubic(unsigned int L) : length_x_(L), length_y_(L), length_z_(L) { init(); }
  cubic(unsigned int Lx, unsigned int Ly, unsigned int Lz) :
    length_x_(Lx), length_y_(Ly), length_z_(Lz) { init(); }
  void init() {
    neighbors_.resize(num_sites());
    source_.resize(num_bonds());
    target_.resize(num_bonds());
    site_phase_.resize(num_sites());
    bond_phase_.resize(num_bonds());
    for (unsigned int s = 0; s < num_sites(); ++s) {
      int x, y, z;
      boost::tie(x, y, z) = index2xyz(s);
      site_phase_[s] = 2.0 * ((x + y + z) % 2) - 1.0;
      for (unsigned int k = 0; k < 6; ++k) {
        int d = 1 - 2 * (k / 3);
        if ((k % 3) == 0) {
          neighbors_[s][k] = xyz2index(x + d, y, z);
        } else if ((k % 3) == 1) {
          neighbors_[s][k] = xyz2index(x, y + d, z);
        } else {
          neighbors_[s][k] = xyz2index(x, y, z + d);
        }
      }
    }
    for (unsigned int b = 0; b < num_bonds(); ++b) {
      unsigned int s = b / 3;
      int x, y, z;
      boost::tie(x, y, z) = index2xyz(s);
      unsigned int t;
      if (b % 3 == 0) {
        t = xyz2index(x + 1, y, z);
        bond_phase_[b] = 2.0 * ((b / 3) % 2) - 1.0;
      } else if (b % 3 == 1) {
        t = xyz2index(x, y + 1, z);
        bond_phase_[b] = 2.0 * ((b / length_x_ / 3) % 2) - 1.0;
      } else {
        t = xyz2index(x, y, z + 1);
        bond_phase_[b] = 2.0 * ((b / length_x_ / length_y_ / 3) % 2) - 1.0;
      }
      source_[b] = s;
      target_[b] = t;
    }
  }
  unsigned int get_length_x() const { return length_x_; }
  unsigned int get_length_y() const { return length_y_; }
  unsigned int get_length_z() const { return length_z_; }
  unsigned int num_sites() const { return length_x_ * length_y_ * length_z_; }
  unsigned int num_bonds() const { return 3 * num_sites(); }
  unsigned int source(unsigned int b) const { return source_[b]; }
  unsigned int target(unsigned int b) const { return target_[b]; }
  unsigned int num_neighbors() const { return 6; }
  unsigned int neighbor(unsigned int s, unsigned int k) const { return neighbors_[s][k]; }
  double site_phase(unsigned int s) const { return site_phase_[s]; }
  double bond_phase(unsigned int b) const { return bond_phase_[b]; }
protected:
  boost::tuple<int, int, int> index2xyz(unsigned int s) const {
    return boost::make_tuple(s % length_x_, (s / length_x_) % length_y_,
                             (s / (length_x_ * length_y_)));
  }
  unsigned int xyz2index(int x, int y, int z) const {
    x += length_x_;
    y += length_y_;
    z += length_z_;
    return x % length_x_ + (y % length_y_) * length_x_ + (z % length_z_) * length_x_ * length_y_;
  }
private:
  unsigned int length_x_, length_y_, length_z_;
  std::vector<boost::array<unsigned int, 6> > neighbors_;
  std::vector<unsigned int> source_;
  std::vector<unsigned int> target_;
  std::vector<double> site_phase_;
  std::vector<double> bond_phase_;
};

} // end namespace lattice

#endif // LATTICE_CUBIC_HPP
