/*****************************************************************************
*
* Copyright (C) 2015-2020 by Synge Todo <wistaria@phy.s.u-tokyo.ac.jp>
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
#include <standards/simpson.hpp>
#include <lattice/graph.hpp>

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
    k_ = 1.0 / (sinh(2 * beta * Jx) * sinh(2 * beta * Jy));
  }
  FVAR operator()(T t) const {
    using std::cos; using std::log; using std::sqrt;
    return c_ * log(chab_ + sqrt(1 + k_ * k_ - 2 * k_ * cos(2 * t)) / k_);
  }
  T c_;
  FVAR chab_, k_;
};

template<typename T, typename FVAR>
functor<T, FVAR> func(T Jx, T Jy, FVAR beta) { return functor<T, FVAR>(Jx, Jy, beta); }
  
}

template<typename T>
inline std::tuple<T, T, T> infinite(T Jx, T Jy, T beta) {
  typedef T real_t;
  if (Jx <= 0 || Jy <= 0)
    throw(std::invalid_argument("Jx and Jy should be positive"));
  if (beta <= 0)
    throw(std::invalid_argument("beta should be positive"));
  real_t pi = boost::math::constants::pi<real_t>();
  auto beta_fvar = boost::math::differentiation::make_fvar<real_t, 2>(beta);
  auto logZ = log(real_t(2)) / 2 + standards::simpson_1d(func(Jx, Jy, beta_fvar), real_t(0), pi, 8);
  real_t d2 = logZ.derivative(2);
  for (unsigned long n = 16; n < 1024; n *= 2) {
    logZ = log(real_t(2)) / 2 + standards::simpson_1d(func(Jx, Jy, beta_fvar), real_t(0), pi, n);
    if (beta * beta * abs((logZ.derivative(2) - d2) / logZ.derivative(0)) <
        2 * std::numeric_limits<real_t>::epsilon()) break;
    d2 = logZ.derivative(2);
  }
  real_t free_energy = - logZ.derivative(0) / beta;
  real_t energy = - logZ.derivative(1);
  real_t specific_heat = beta * beta * logZ.derivative(2);
  return std::make_tuple(free_energy, energy, specific_heat);
}

template<typename I, typename T>
inline std::tuple<T, T, T> finite(I Lx, I Ly, T Jx, T Jy, T beta) {
  typedef I int_t;
  typedef T real_t;
  if (Lx <= 0 || Ly <= 0)
    throw(std::invalid_argument("Lx and Ly should be positive"));
  if (Jx <= 0 || Jy <= 0)
    throw(std::invalid_argument("Jx and Jy should be positive"));
  if (beta <= 0)
    throw(std::invalid_argument("beta should be positive"));
  real_t pi = boost::math::constants::pi<real_t>();
  auto beta_fvar = boost::math::differentiation::make_fvar<real_t, 2>(beta);
  auto a = beta_fvar * Jx;
  auto b = beta_fvar * Jy;
  decltype(beta_fvar) lp0(0), lp1(0), lp2(0), lp3(0);
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
  if (sinh(2*a) * sinh(2*b) < real_t(1)) {
    logZ += (lp0 + log(1 + exp(lp1-lp0) + exp(lp2-lp0) - exp(lp3-lp0))) / (Lx * Ly);
  } else {
    logZ += (lp0 + log(1 + exp(lp1-lp0) + exp(lp2-lp0) + exp(lp3-lp0))) / (Lx * Ly);
  }
  real_t free_energy = - logZ.derivative(0) / beta;
  real_t energy = - logZ.derivative(1);
  real_t specific_heat = beta * beta * logZ.derivative(2);
  return std::make_tuple(free_energy, energy, specific_heat);
}

template<typename I, typename T>
inline std::tuple<T, T, T> finite_count(I Lx, I Ly, T Jx, T Jy, T beta) {
  typedef I int_t;
  typedef T real_t;
  if (Lx <= 0 || Ly <= 0)
    throw(std::invalid_argument("Lx and Ly should be positive"));
  if (Jx <= 0 || Jy <= 0)
    throw(std::invalid_argument("Jx and Jy should be positive"));
  if (beta <= 0)
    throw(std::invalid_argument("beta should be positive"));
  auto beta_fvar = boost::math::differentiation::make_fvar<real_t, 2>(beta);

  auto basis = lattice::basis::simple(2);
  auto unitcell = lattice::unitcell(2);
  unitcell.add_site(lattice::coordinate(0, 0), 0);
  unitcell.add_bond(0, 0, lattice::offset(1, 0), 0);
  unitcell.add_bond(0, 0, lattice::offset(0, 1), 1);
  auto graph = lattice::graph(basis, unitcell, lattice::extent(Lx, Ly));
  
  std::vector<real_t> J(2);
  J[0] = Jx;
  J[1] = Jy;
  real_t gs_energy = -(Jx + Jy) * (Lx * Ly);
  int_t num_states = 1 << graph.num_sites();
  decltype(beta_fvar) Z(0);
  for (int_t  c = 0; c < num_states; ++c) {
    real_t energy = 0;
    for (int_t b = 0; b < graph.num_bonds(); ++b) {
      real_t ci = real_t(2) * ((c >> graph.source(b)) & 1) - 1;
      real_t cj = real_t(2) * ((c >> graph.target(b)) & 1) - 1;
      energy -= J[graph.bond_type(b)] * ci * cj;
    }
    Z += exp(- beta_fvar * (energy - gs_energy));
  }
  auto logZ = (log(Z) - beta_fvar * gs_energy) / (Lx * Ly);
  real_t free_energy = - logZ.derivative(0) / beta;
  real_t energy = - logZ.derivative(1);
  real_t specific_heat = beta * beta * logZ.derivative(2);
  return std::make_tuple(free_energy, energy, specific_heat);
}

} // end namespace square
} // end namespace free_energy
} // end namespace ising
