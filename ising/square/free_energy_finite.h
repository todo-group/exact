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

#ifndef ISING_SQUARE_FREE_ENERGY_FINITE_H
#define ISING_SQUARE_FREE_ENERGY_FINITE_H

#include <log_sum_exp/exp_number.hpp>
#include <vector>

namespace {
  
inline log_sum_exp::exp_double cosh_value(double x) {
  return (log_sum_exp::exp_value(x) + log_sum_exp::exp_value(-x)) / 2;
}

inline log_sum_exp::exp_double sinh_value(double x) {
  return (log_sum_exp::exp_value(x) - log_sum_exp::exp_value(-x)) / 2;
}

}

namespace ising {
namespace square {

inline log_sum_exp::exp_double partition_function(double beta, double Jx, double Jy, int Lx, int Ly) {
  double a = beta * Jx;
  double b = beta * Jy;
  std::vector<double> gamma(2 * Lx);
  for (int k = 0; k < 2 * Lx; ++k) {
    log_sum_exp::exp_double cosh_g =
      (cosh_value(2*a) * cosh_value(2*b) - cos(M_PI*k/Lx) * sinh_value(2*b)) / sinh_value(2*a);
    gamma[k] = log(cosh_g + sqrt(cosh_g * cosh_g - 1));
  }
  if (sinh_value(2*a) * sinh_value(2*b) > 1) gamma[0] = -gamma[0];
  log_sum_exp::exp_double p0(1), p1(1), p2(1), p3(1);
  for (int k = 1; k <= Lx; ++k) {
    p0 *= 2 * cosh_value(Ly * gamma[2*k-1] / 2);
    p1 *= 2 * sinh_value(Ly * gamma[2*k-1] / 2);
    p2 *= 2 * cosh_value(Ly * gamma[2*k-2] / 2);
    p3 *= 2 * sinh_value(Ly * gamma[2*k-2] / 2);
  }
  log_sum_exp::exp_double z = 0.5 * pow(2 * sinh_value(2*a), Lx*Ly/2) * (p0 + p1 + p2 - p3);
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

#endif // ISING_SQUARE_FREE_ENERGY_FINITE_H
