/*****************************************************************************
*
* Copyright (C) 2015-2016 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Calculating free energy density of triangular lattice Ising model in the thermodynamic limit

// reference: R. M. F. Houtappel, Physica 16, 425 (1950)

#ifndef ISING_TRIANGLE_FREE_ENERGY_H
#define ISING_TRIANGLE_FREE_ENERGY_H

#include <boost/throw_exception.hpp>
#include <cmath>
#include <stdexcept>

namespace {

// ch2a = cosh(2 * beta * Ja), sh2a = sinh(2 * beta * Ja), etc
inline double f(double omega1, double omega2, double ch2a, double ch2b, double ch2c,
                double sh2a, double sh2b, double sh2c) {
  return std::log(ch2a * ch2b * ch2c + sh2a * sh2b * sh2c
                  - sh2a * std::cos(omega1) - sh2b * std::cos(omega2)
                  - sh2c * std::cos(omega1 + omega2));
}

}

namespace ising {
namespace triangle {
  
inline double free_energy_density(double beta, double Ja, double Jb, double Jc, int Nint = 100) {
  double ch2a = std::cosh(2 * beta * Ja);
  double ch2b = std::cosh(2 * beta * Jb);
  double ch2c = std::cosh(2 * beta * Jc);
  double sh2a = std::sinh(2 * beta * Ja);
  double sh2b = std::sinh(2 * beta * Jb);
  double sh2c = std::sinh(2 * beta * Jc);
  double domega = 2 * M_PI / Nint;
  double g = 0.0;
  // i == 0
  {
    double omega1 = 0;
    //   j == 0 || i == Nint
    g += f(omega1, 0, ch2a, ch2b, ch2c, sh2a, sh2b, sh2c)
      + f(omega1, 2 * M_PI, ch2a, ch2b, ch2c, sh2a, sh2b, sh2c);
    //   j = 1/2 ... Nint-1/2
    for (int j = 0; j < Nint; ++j) {
      double omega2 = domega * (j + 0.5);
      g += 4 * f(omega1, omega2, ch2a, ch2b, ch2c, sh2a, sh2b, sh2c);
    }
    //   j = 1 ... Nint-1
    for (int j = 1; j < Nint; ++j) {
      double omega2 = domega * j;
      g += 2 * f(omega1, omega2, ch2a, ch2b, ch2c, sh2a, sh2b, sh2c);
    }
  }
  // i = 1/2 ... Nint-1/2
  for (int i = 0; i < Nint; ++i) {
    double omega1 = domega * (i + 0.5);
    //   j == 0 || i == Nint
    g += 4 * f(omega1, 0, ch2a, ch2b, ch2c, sh2a, sh2b, sh2c)
      + 4 * f(omega1, 2 * M_PI, ch2a, ch2b, ch2c, sh2a, sh2b, sh2c);
    //   j = 1/2 ... Nint-1/2
    for (int j = 0; j < Nint; ++j) {
      double omega2 = domega * (j + 0.5);
      g += 16 * f(omega1, omega2, ch2a, ch2b, ch2c, sh2a, sh2b, sh2c);
    }
    //   j = 1 ... Nint-1
    for (int j = 1; j < Nint; ++j) {
      double omega2 = domega * j;
      g += 8 * f(omega1, omega2, ch2a, ch2b, ch2c, sh2a, sh2b, sh2c);
    }
  }
  // i = 1 ... Nint-1
  for (int i = 1; i < Nint; ++i) {
    double omega1 = domega * i;
    //   j == 0 || i == Nint
    g += 2 * f(omega1, 0, ch2a, ch2b, ch2c, sh2a, sh2b, sh2c)
      + 2 * f(omega1, 2 * M_PI, ch2a, ch2b, ch2c, sh2a, sh2b, sh2c);
    //   j = 1/2 ... Nint-1/2
    for (int j = 0; j < Nint; ++j) {
      double omega2 = domega * (j + 0.5);
      g += 8 * f(omega1, omega2, ch2a, ch2b, ch2c, sh2a, sh2b, sh2c);
    }
    //   j = 1 ... Nint-1
    for (int j = 1; j < Nint; ++j) {
      double omega2 = domega * j;
      g += 4 * f(omega1, omega2, ch2a, ch2b, ch2c, sh2a, sh2b, sh2c);
    }
  }
  // i == Nint
  {
    double omega1 = 2 * M_PI;
    //   j == 0 || i == Nint
    g += f(omega1, 0, ch2a, ch2b, ch2c, sh2a, sh2b, sh2c)
      + f(omega1, 2 * M_PI, ch2a, ch2b, ch2c, sh2a, sh2b, sh2c);
    //   j = 1/2 ... Nint-1/2
    for (int j = 0; j < Nint; ++j) {
      double omega2 = domega * (j + 0.5);
      g += 4 * f(omega1, omega2, ch2a, ch2b, ch2c, sh2a, sh2b, sh2c);
    }
    //   j = 1 ... Nint-1
    for (int j = 1; j < Nint; ++j) {
      double omega2 = domega * j;
      g += 2 * f(omega1, omega2, ch2a, ch2b, ch2c, sh2a, sh2b, sh2c);
    }
  }

  g *= domega * domega / 36;
  g /= 8 * M_PI * M_PI;
  return - (std::log(2.0) + g) / beta;
}

} // end namespace triangle
} // end namespace ising

#endif // ISING_TRIANGLE_FREE_ENERGY_H
