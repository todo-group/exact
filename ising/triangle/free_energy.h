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

#include "simpson_integration.h"
#include <boost/throw_exception.hpp>
#include <cmath>
#include <stdexcept>

namespace {

struct func {
  func(double beta, double Ja, double Jb, double Jc) {
    c_ = 1.0 / (8 * M_PI * M_PI);
    sh2a_ = std::sinh(2 * beta * Ja);
    sh2b_ = std::sinh(2 * beta * Jb);
    sh2c_ = std::sinh(2 * beta * Jc);
    cshabc_ = std::cosh(2 * beta * Ja) * std::cosh(2 * beta * Jb) * std::cosh(2 * beta * Jc)
      + sh2a_ * sh2b_ * sh2c_;
  }
  double operator()(double t1, double t2) const {
    return c_ * std::log(cshabc_ - sh2a_ * std::cos(t1) - sh2b_ * std::cos(t2)
                         - sh2c_ * std::cos(t1 + t2));
  }
  double c_, sh2a_, sh2b_, sh2c_, cshabc_;
};

}

namespace ising {
namespace triangle {
  
inline double free_energy_density(double beta, double Ja, double Jb, double Jc, int Nint) {
  if (beta <= 0)
    boost::throw_exception(std::invalid_argument("beta should be positive"));
  func f(beta, Ja, Jb, Jc);
  return - (std::log(2.0) + simpson_integration_2d(f, 0, 0, 2 * M_PI, 2 * M_PI, Nint, Nint)) / beta;
}

} // end namespace triangle
} // end namespace ising

#endif // ISING_TRIANGLE_FREE_ENERGY_H
