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

#ifndef ISING_SQUARE_FREE_ENERGY_H
#define ISING_SQUARE_FREE_ENERGY_H

#include <boost/throw_exception.hpp>
#include <cmath>
#include <stdexcept>

namespace {

// ch2a = cosh(2 * beta * Ja), sh2a = sinh(2 * beta * Ja), etc
inline double f(double theta, double ch2a, double ch2b, double sh2a, double sh2b) {
  double k = 1.0 / (sh2a * sh2b);
  return std::log(ch2a * ch2b + std::sqrt(1 + k * k - 2 * k * std::cos(2 * theta)) / k);
}

}

namespace ising {
namespace square {
  
inline double free_energy_density(double beta, double Ja, double Jb, unsigned int Nint) {
  if (beta <= 0)
    boost::throw_exception(std::invalid_argument("beta should be positive"));
  if (Nint == 0 || (Nint % 2) != 0)
    boost::throw_exception(std::invalid_argument("Nint should be positive and a multiple of two"));
  double ch2a = std::cosh(2 * beta * Ja);
  double ch2b = std::cosh(2 * beta * Jb);
  double sh2a = std::sinh(2 * beta * Ja);
  double sh2b = std::sinh(2 * beta * Jb);
  double dt = M_PI / Nint;
  double g = 0.0;
  // i == 0 || i == Nint
  g += f(0, ch2a, ch2b, sh2a, sh2b) + f(M_PI, ch2a, ch2b, sh2a, sh2b);
  // i = 1/2 ... Nint-1/2
  for (int i = 0; i < Nint; ++i) {
    double t = dt * (i + 0.5);
    g += 4 * f(t, ch2a, ch2b, sh2a, sh2b);
  }
  // i = 1 ... Nint-1
  for (int i = 1; i < Nint; ++i) {
    double t = dt * i;
    g += 2 * f(t, ch2a, ch2b, sh2a, sh2b);
  }
  g *= dt / 6;
  g /= 2 * M_PI;
  return - (std::log(2.0) / 2 + g) / beta;
}

} // end namespace square
} // end namespace ising

#endif // ISING_SQUARE_FREE_ENERGY_H
