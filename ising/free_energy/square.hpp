/*****************************************************************************
*
* Copyright (C) 2015-2021 by Synge Todo <wistaria@phy.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Free energy, energy, and specific heat of square lattice Ising model

#pragma once

#include <cmath>
#include <stdexcept>
#include <boost/math/constants/constants.hpp>
#include <boost/math/differentiation/autodiff.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <standards/simpson.hpp>
#include <lattice/graph.hpp>
#include "common.hpp"

namespace ising {
namespace free_energy {
namespace square {

namespace {

template<typename T, typename FVAR>
struct functor {
  functor(T Jx, T Jy, FVAR beta) {
    using std::cosh; using std::sinh;
    c_ = 1 / (2 * boost::math::constants::pi<T>());
    chab_ = cosh(2 * beta * Jx) * cosh(2 * beta * Jy);
    k_ = 1 / (sinh(2 * beta * Jx) * sinh(2 * beta * Jy));
    std::cerr << k_ << ' ' << chab_ << ' ' << std::endl;
  }
  FVAR operator()(T t) const {
    using std::cos; using std::log; using std::sqrt;
    // auto x = k_ - 1;
    return c_ * log(chab_ + sqrt(1 + k_ * k_ - 2 * k_ * cos(2 * t)) / k_);
    // return c_ * log(chab_ + sqrt(x * x + 4 * (1+x) * sin(t) * sin(t)) / k_);
  }
  T c_;
  FVAR chab_, k_;
};

template<typename T, typename FVAR>
functor<T, FVAR> func(T Jx, T Jy, FVAR beta) { return functor<T, FVAR>(Jx, Jy, beta); }
  
}

template<typename T, typename U>
inline U infinite(T Jx, T Jy, U beta) {
  const unsigned long max_n = 1 << 16;
  typedef T real_t;
  if (Jx <= 0 || Jy <= 0)
    throw(std::invalid_argument("Jx and Jy should be positive"));
  if (beta <= 0)
    throw(std::invalid_argument("beta should be positive"));
  real_t pi = boost::math::constants::pi<real_t>();
  // auto logZ = log(real_t(2)) / 2 + 2 * standards::simpson_1d(func(Jx, Jy, beta), real_t(0), pi / 2, 8);
  boost::math::quadrature::tanh_sinh<real_t> integrator;
  auto logZ = log(real_t(2)) / 2 + 2 * integrator.integrate(func(Jx, Jy, beta), 0, pi / 2);
  // auto pz = logZ;
  // auto pe = abs(beta * beta * logZ.derivative(2) / logZ);
  // for (unsigned long n = 16; n <= max_n; n *= 2) {
  //   logZ = log(real_t(2)) / 2 + 2 * standards::simpson_1d(func(Jx, Jy, beta), real_t(0), pi / 2, n);
  //   auto pn = abs(beta * beta * (logZ - pz).derivative(2) / logZ);
  //   if (pn < 2 * std::numeric_limits<real_t>::epsilon() || pn > pe) break;
  //   if (n == max_n) std::cerr << "Warning: integration not converge with n = " << max_n << std::endl;
  //   pz = logZ;
  //   pe = pn;
  // }
  return - logZ / beta;
}

template<typename I, typename T, typename U>
inline U finite(I Lx, I Ly, T Jx, T Jy, U beta) {
  typedef I int_t;
  typedef T real_t;
  typedef U value_t;
  if (Lx <= 0 || Ly <= 0)
    throw(std::invalid_argument("Lx and Ly should be positive"));
  if (Jx <= 0 || Jy <= 0)
    throw(std::invalid_argument("Jx and Jy should be positive"));
  real_t pi = boost::math::constants::pi<real_t>();
  auto a = beta * Jx;
  auto b = beta * Jy;
  value_t lp0(0), lp1(0), lp2(0), lp3(0);
  for (int_t k = 0; k < 2 * Ly; ++k) {
    auto cosh_g = (cosh(2*a) * cosh(2*b) - cos(pi * k / Ly) * sinh(2*b)) / sinh(2*a);
    auto gamma_abs = log(cosh_g + sqrt(cosh_g * cosh_g - 1));
    if ((k & 1) == 1) {
      lp0 += (Lx * gamma_abs/2) + log(1 + exp(-(Lx * gamma_abs)));
      lp1 += (Lx * gamma_abs/2) + log(1 - exp(-(Lx * gamma_abs)));
    } else {
      lp2 += (Lx * gamma_abs/2) + log(1 + exp(-(Lx * gamma_abs)));
      lp3 += (Lx * gamma_abs/2) + log(1 - exp(-(Lx * gamma_abs)));
    }
  }
  auto logZ = -log(real_t(2)) / (Lx * Ly) + a + real_t(1)/2 * log(1 - exp(-4*a));
  if (sinh(2*a) * sinh(2*b) < value_t(1)) {
    logZ += (lp0 + log(1 + exp(lp1-lp0) + exp(lp2-lp0) - exp(lp3-lp0))) / (Lx * Ly);
  } else {
    logZ += (lp0 + log(1 + exp(lp1-lp0) + exp(lp2-lp0) + exp(lp3-lp0))) / (Lx * Ly);
  }
  return - logZ / beta;
}

template<typename T, typename I>
inline T finite_tc(I Lx, I Ly) {
  typedef I int_t;
  typedef T real_t;
  if (Lx <= 0 || Ly <= 0)
    throw(std::invalid_argument("Lx and Ly should be positive"));
  real_t pi = boost::math::constants::pi<real_t>();
  real_t beta = log(sqrt(real_t(2)) + 1)/2;
  real_t lp0(0), lp1(0);
  for (int_t k = 1; k < 2 * Ly; k += 2) {
    real_t c = (2 - cos(pi * k / Ly));
    real_t cosh_g = c;
    real_t gamma_abs = log(cosh_g + sqrt(cosh_g * cosh_g - 1));
    lp0 += (Lx * gamma_abs/2) + log(1 + exp(-(Lx * gamma_abs)));
    lp1 += (Lx * gamma_abs/2) + log(1 - exp(-(Lx * gamma_abs)));
  }
  real_t logZ = -log(real_t(2)) / (Lx * Ly) + beta + real_t(1)/2 * log(1 - exp(-4*beta))
    + (lp0 + 2 * lp1) / (Lx * Ly)
    + (lp0 + log(1 + 2 * exp(lp1-lp0))) / (Lx * Ly);
  return - logZ / beta;
}
// template<typename T, typename I>
// inline std::tuple<T, T, T> finite_tc(I Lx, I Ly) {
//   typedef I int_t;
//   typedef T real_t;
//   if (Lx <= 0 || Ly <= 0)
//     throw(std::invalid_argument("Lx and Ly should be positive"));
//   auto pi = boost::math::constants::pi<real_t>();
//   auto beta = log(sqrt(real_t(2))+1)/2;
//   auto beta_fvar = boost::math::differentiation::make_fvar<real_t, 2>(beta);
//   decltype(beta_fvar) lp0(0), lp1(0);
//   auto x = boost::math::differentiation::make_fvar<real_t, 2>(0);
//   auto one = real_t(1);
//   auto two = real_t(2);
//   auto cosh_a = sqrt(2) + 2 * x + 2 * sqrt(2) * x * x;
//   auto sinh_a = 1 + 2 * sqrt(2) * x + 2 * x * x;
//   std::cout << cosh_a << ' ' << sinh_a << std::endl;
//   for (int_t k = 1; k < 2 * Ly; k += 2) {
//     auto c = (2 - cos(pi * k / Ly));
//     auto cosh_g = c + 8 * x * x;
//     auto gamma_abs = log(cosh_g + sqrt(cosh_g * cosh_g - 1));
//     std::cout << k << ' ' << cosh_g << ' ' << sqrt(cosh_g * cosh_g - 1) << ' ' << gamma_abs << ' ' << log(1 + exp(-(Lx * gamma_abs))) << ' ' << log(1 - exp(-(Lx * gamma_abs))) << std::endl;
//     lp0 += (Lx * gamma_abs/2) + log(1 + exp(-(Lx * gamma_abs)));
//     lp1 += (Lx * gamma_abs/2) + log(1 - exp(-(Lx * gamma_abs)));
//   }
//   std::cout << lp0 << ' ' << lp1 << std::endl;
//   auto logZ = -log(real_t(2)) / (Lx * Ly) + beta_fvar + real_t(1)/2 * log(1 - exp(-4*beta_fvar))
//     + (lp0 + log(1 + 2 * exp(lp1-lp0))) / (Lx * Ly);
//   std::cout << logZ << std::endl;
//   real_t free_energy = - logZ.derivative(0) / beta;
//   real_t energy = - logZ.derivative(1);
//   real_t specific_heat = beta * beta * logZ.derivative(2);
//   return std::make_tuple(free_energy, energy, specific_heat);
// }

template<typename I, typename T, typename U>
inline U finite_count(I Lx, I Ly, T Jx, T Jy, U beta) {
  typedef I int_t;
  typedef T real_t;
  typedef U value_t;
  if (Lx <= 0 || Ly <= 0)
    throw(std::invalid_argument("Lx and Ly should be positive"));
  if (Jx <= 0 || Jy <= 0)
    throw(std::invalid_argument("Jx and Jy should be positive"));

  auto basis = lattice::basis::simple(2);
  auto unitcell = lattice::unitcell(2);
  unitcell.add_site(lattice::coordinate(0, 0), 0);
  unitcell.add_bond(0, 0, lattice::offset(1, 0), 0);
  unitcell.add_bond(0, 0, lattice::offset(0, 1), 1);
  auto graph = lattice::graph(basis, unitcell, lattice::extent(Lx, Ly));
  if (graph.num_sites() >= 32) throw(std::invalid_argument("too large system size"));
  
  std::vector<real_t> J(2);
  J[0] = Jx;
  J[1] = Jy;
  real_t gs_energy = -(Jx + Jy) * (Lx * Ly);
  int_t num_states = 1 << graph.num_sites();
  value_t Z(0);
  for (int_t c = 0; c < num_states; ++c) {
    real_t energy = 0;
    for (int_t b = 0; b < graph.num_bonds(); ++b) {
      real_t ci = real_t(2) * ((c >> graph.source(b)) & 1) - 1;
      real_t cj = real_t(2) * ((c >> graph.target(b)) & 1) - 1;
      energy -= J[graph.bond_type(b)] * ci * cj;
    }
    Z += exp((gs_energy - energy) * beta);
  }
  return (-log(Z)/beta + gs_energy) / (Lx * Ly);
}

} // end namespace square
} // end namespace free_energy
} // end namespace ising
