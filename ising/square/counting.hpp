/*****************************************************************************
*
* Copyright (C) 2015-2017 by Synge Todo <wistaria@phy.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Calculating free energy density of square lattice Ising model by exact counting

#include <lattice/square.hpp>
#include <lse/exp_number.hpp>
#include <boost/throw_exception.hpp>
#include <stdexcept>
#include <vector>

#ifndef ISING_SQUARE_COUNTING_HPP
#define ISING_SQUARE_COUNTING_HPP

namespace ising {
namespace square {

struct counting {
public:
  template<typename LATTICE>
  static double free_energy(double beta, LATTICE const& lat, std::vector<double> const& inter,
                            std::vector<double> const& field = std::vector<double>(0)) {
    if (beta <= 0)
      boost::throw_exception(std::invalid_argument("beta should be positive"));
    if (lat.num_sites() > 30)
      boost::throw_exception(std::invalid_argument("too large lattice"));
    if (inter.size() != lat.num_bonds())
      boost::throw_exception(std::invalid_argument("inconsitent table size of interaction"));
    if (field.size() > 0 && field.size() != lat.num_sites())
      boost::throw_exception(std::invalid_argument("inconsitent table size of external field"));
    unsigned long num_states = 1 << lat.num_sites();
    lse::exp_double sum = 0;
    for (unsigned long c = 0; c < num_states; ++c) {
      lse::exp_double weight = 1;
      for (unsigned int b = 0; b < lat.num_bonds(); ++b) {
        int ci = (c >> lat.source(b)) & 1;
        int cj = (c >> lat.target(b)) & 1;
        weight *= lse::exp_value(beta * inter[b] * (1 - 2 * (ci ^ cj)));
      }
      if (field.size()) {
        for (unsigned int s = 0; s < lat.num_sites(); ++s) {
          int c0 = (c >> s) & 1;
          weight *= lse::exp_value(beta * field[s] * (1 - 2 * c0));
        }
      }
      sum += weight;
    }
    return -log(sum) / beta;
  }

  template<typename LATTICE>
  static double free_energy(double beta, LATTICE const& lat, double J, double H = 0.0) {
    std::vector<double> inter(lat.num_bonds(), J);
    std::vector<double> field((H != 0.0) ? lat.num_sites() : 0, H);
    return free_energy(beta, lat, inter, field);
  }

  template<typename LATTICE>
  static double free_energy_density(double beta, LATTICE const& lat,
                                    std::vector<double> const& inter,
                                    std::vector<double> const& field = std::vector<double>(0)) {
    return free_energy(beta, lat, inter, field) / lat.num_sites();
  }

  template<typename LATTICE>
  static double free_energy_density(double beta, LATTICE const& lat, double J, double H = 0.0) {
    return free_energy(beta, lat, J, H) / lat.num_sites();
  }
};

} // end namespace square
} // end namespace ising

#endif // ISING_SQUARE_COUNTING_HPP
