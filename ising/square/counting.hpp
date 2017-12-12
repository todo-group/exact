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
  static double free_energy(double beta, lattice::square const& lat,
                            std::vector<double> interaction) {
    if (beta <= 0)
      boost::throw_exception(std::invalid_argument("beta should be positive"));
    if (lat.num_sites() > 30)
      boost::throw_exception(std::invalid_argument("too large lattice"));
    unsigned long num_states = 1 << lat.num_sites();
    lse::exp_double sum = 0;
    for (unsigned long s = 0; s < num_states; ++s) {
      lse::exp_double weight = 1;
      for (unsigned int b = 0; b < lat.num_bonds(); ++b) {
        int si = (s >> lat.source(b)) & 1;
        int sj = (s >> lat.target(b)) & 1;
        weight *= lse::exp_value(beta * interaction[b] * (1 - 2 * (si ^ sj)));
      }
      sum += weight;
    }
    return -log(sum) / beta;
  }

  static double free_energy_density(double beta, lattice::square const& lat,
                                    std::vector<double> interaction) {
    return free_energy(beta, lat, interaction) / lat.num_sites();
  }
};

} // end namespace square
} // end namespace ising

#endif // ISING_SQUARE_COUNTING_HPP
