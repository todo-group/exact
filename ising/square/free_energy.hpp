/*****************************************************************************
*
* Copyright (C) 2015-2016 by Synge Todo <wistaria@phy.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Calculating free energy density of square lattice Ising model in the thermodynamic limit

// reference: https://en.wikipedia.org/wiki/Square-lattice_Ising_model

#ifndef ISING_SQUARE_FREE_ENERGY_HPP
#define ISING_SQUARE_FREE_ENERGY_HPP

#include <cmath>
#include <stdexcept>
#include <boost/throw_exception.hpp>
#include <integral/simpson.hpp>

namespace {

struct func {
  func(double beta, double Ja, double Jb) {
    c_ = 1.0 / (2 * M_PI);
    chab_ = std::cosh(2 * beta * Ja) * std::cosh(2 * beta * Jb);
    k_ = 1.0 / (std::sinh(2 * beta * Ja) * std::sinh(2 * beta * Jb));
  }
  double operator()(double t) const {
    return c_ * std::log(chab_ + std::sqrt(1 + k_ * k_ - 2 * k_ * std::cos(2 * t)) / k_);
  }
  double c_, chab_, k_;
};

}

namespace ising {
namespace square {
  
inline double free_energy_density(double beta, double Ja, double Jb, unsigned int Nint) {
  if (beta <= 0)
    boost::throw_exception(std::invalid_argument("beta should be positive"));
  func f(beta, Ja, Jb);
  return - (std::log(2.0) / 2 + integral::simpson_1d(f, 0, M_PI, Nint)) / beta;
}

} // end namespace square
} // end namespace ising

#endif // ISING_SQUARE_FREE_ENERGY_HPP
