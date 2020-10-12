/*****************************************************************************
*
* Copyright (C) 2011-2016 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Calculating free energy density of square lattice Ising model

// reference: B. Kastening, Phys. Rev. E 64, 066106 (2001), wrn:2011/02/10

#ifndef ISING_SQUARE_FREE_ENERGY_FINITE_HPP
#define ISING_SQUARE_FREE_ENERGY_FINITE_HPP

#include <vector>
#include <lse/exp_number.hpp>

namespace {
  
inline lse::exp_double cosh_value(double x) {
  return (lse::exp_value(x) + lse::exp_value(-x)) / 2;
}

inline lse::exp_double sinh_value(double x) {
  return (lse::exp_value(x) - lse::exp_value(-x)) / 2;
}

}

namespace ising {
namespace square {

inline double partition_function(double beta, double Jx, double Jy, int Lx, int Ly) {
  typedef double real_t;
  auto a = beta * Jx;
  auto b = beta * Jy;
  auto a_bar = log((1 + cosh(2*a)) / sinh(2*a)) / 2;
  double p0(1), p1(1), p2(1), p3(1);
  double lp0(0), lp1(0), lp2(0), lp3(0);
  for (int k = 0; k < 2 * Lx; ++k) {
    auto cosh_g =
      (cosh(2*a) * cosh(2*b) - cos(M_PI*k/Lx) * sinh(2*b)) / sinh(2*a);
    auto gamma = log(cosh_g + sqrt(cosh_g * cosh_g - 1));
    // if (k == 0) gamma = 2 * (a_bar - b);
    if ((k & 1) == 0) {
      //p2 *= exp(Ly * gamma / 2) + exp(- (Ly * gamma / 2));
      //p3 *= exp(Ly * gamma / 2) - exp(- (Ly * gamma / 2));
      p2 *= 2 * cosh(Ly * gamma / 2);
      p3 *= 2 * sinh(Ly * gamma / 2);
      lp2 += log(2 * cosh(Ly * gamma / 2));
      lp3 += log(2 * sinh(Ly * gamma / 2));
    } else {
      //p0 *= exp(Ly * gamma / 2) + exp(- (Ly * gamma / 2));
      //p1 *= exp(Ly * gamma / 2) - exp(- (Ly * gamma / 2));
      p0 *= 2 * cosh(Ly * gamma / 2);
      p1 *= 2 * sinh(Ly * gamma / 2);
      lp0 += log(2 * cosh(Ly * gamma / 2));
      lp1 += log(2 * sinh(Ly * gamma / 2));
    }
  }
  auto logZ = - log(real_t(2)) / (Lx * Ly) + a + 0.5 * log(1 - exp(-4*a));
  if (a_bar > b) {
    logZ += (lp0 + log(1 + exp(lp1-lp0) + exp(lp2-lp0) - exp(lp3-lp0))) / (Lx * Ly);
  } else {
    logZ += (lp0 + log(1 + exp(lp1-lp0) + exp(lp2-lp0) + exp(lp3-lp0))) / (Lx * Ly);
  }
  auto z = exp(Lx * Ly * logZ);
  return z;
}

inline double free_energy(double beta, double Jx, double Jy, int Lx, int Ly) {
  return -log(partition_function(beta, Jx, Jy, Lx, Ly)) / beta;
}

inline double free_energy_density(double beta, double Jx, double Jy, int Lx, int Ly) {
  return free_energy(beta, Jx, Jy, Lx, Ly) / (Lx * Ly);
}

} // end namespace square
} // end namespace ising

#endif // ISING_SQUARE_FREE_ENERGY_FINITE_HPP
