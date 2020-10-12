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

#ifndef ISING_SQUARE_FINITE_HPP
#define ISING_SQUARE_FINITE_HPP

#include <vector>
#include <cmath>
#include <boost/math/differentiation/autodiff.hpp>
// #include <lse/exp_number.hpp>

// namespace {
  
// inline lse::exp_double cosh_value(double x) {
//   return (lse::exp_value(x) + lse::exp_value(-x)) / 2;
// }

// inline lse::exp_double sinh_value(double x) {
//   return (lse::exp_value(x) - lse::exp_value(-x)) / 2;
// }

// }

namespace ising {
namespace square {

template <typename FVAR>
inline double partition_function_impl(const FVAR& beta, double Jx, double Jy, int Lx, int Ly) {
  auto a = beta * Jx;
  auto b = beta * Jy;
  std::vector<FVAR> gamma(2 * Lx);
  for (int k = 0; k < 2 * Lx; ++k) {
    auto cosh_g =
      (cosh(2*a) * cosh(2*b) - cos(M_PI*k/Lx) * sinh(2*b)) / sinh(2*a);
    gamma[k] = log(cosh_g + sqrt(cosh_g * cosh_g - 1));
  }
  if (sinh(2*a) * sinh(2*b) > 1) gamma[0] = -gamma[0];
  FVAR p0(1), p1(1), p2(1), p3(1);
  for (int k = 1; k <= Lx; ++k) {
    p0 *= 2 * cosh(Ly * gamma[2*k-1] / 2);
    p1 *= 2 * sinh(Ly * gamma[2*k-1] / 2);
    p2 *= 2 * cosh(Ly * gamma[2*k-2] / 2);
    p3 *= 2 * sinh(Ly * gamma[2*k-2] / 2);
  }
  auto z = 0.5 * pow(2 * sinh(2*a), Lx*Ly/2) * (p0 + p1 + p2 - p3);
  auto f = -log(z) / beta / (Lx * Ly);
  // auto e = 
  std::cout << f.derivative(0) << ' ' << f.derivative(1) << ' ' << f.derivative(2) << std::endl;
  return z.derivative(0);
}

inline double partition_function(double beta_in, double Jx, double Jy, int Lx, int Ly) {
  using namespace boost::math::differentiation;
  constexpr unsigned Order = 2;
  auto const beta = make_fvar<double, Order>(beta_in);
  return partition_function_impl(beta, Jx, Jy, Lx, Ly);
  
}
  
inline double free_energy(double beta, double Jx, double Jy, int Lx, int Ly) {
  return -log(partition_function(beta, Jx, Jy, Lx, Ly)) / beta;
}

inline double free_energy_density(double beta, double Jx, double Jy, int Lx, int Ly) {
  return free_energy(beta, Jx, Jy, Lx, Ly) / (Lx * Ly);
}

} // end namespace square
} // end namespace ising

#endif // ISING_SQUARE_FINITE_HPP
