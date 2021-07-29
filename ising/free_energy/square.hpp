/*
   Copyright (C) 2015-2021 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

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
  typedef T real_t;
  typedef FVAR fvar_t;
  functor(real_t Jx, real_t Jy, fvar_t beta) {
    using std::cosh; using std::sinh;
    c_ = 1 / (2 * boost::math::constants::pi<real_t>());
    chab_ = cosh(2 * beta * Jx) * cosh(2 * beta * Jy);
    k_ = 1 / (sinh(2 * beta * Jx) * sinh(2 * beta * Jy));
  }
  fvar_t operator()(real_t t) const {
    using std::cos; using std::log; using std::sqrt;
    return c_ * log(chab_ + sqrt(1 + k_ * k_ - 2 * k_ * cos(2 * t)) / k_);
  }
  real_t c_;
  fvar_t chab_, k_;
};

template<typename T, std::size_t Order>
struct functor<T, boost::math::differentiation::detail::fvar<T, Order>> {
  typedef T real_t;
  typedef boost::math::differentiation::detail::fvar<real_t, Order> fvar_t;
  functor(real_t Jx, real_t Jy, fvar_t beta) {
    using std::cosh; using std::sinh;
    c_ = 1 / (2 * boost::math::constants::pi<real_t>());
    chab_ = cosh(2 * beta * Jx) * cosh(2 * beta * Jy);
    k_ = 1 / (sinh(2 * beta * Jx) * sinh(2 * beta * Jy));
  }
  fvar_t operator()(real_t t) const {
    using std::abs; using std::cos; using std::log; using std::sqrt;
    auto r = sqrt(1 + k_ * k_ - 2 * k_ * cos(2 * t)) / k_;
    if (abs(r.derivative(1)) > std::numeric_limits<real_t>::max()) {
      auto x = boost::math::differentiation::make_fvar<real_t, 2>(real_t(0));
      r = fvar_t(r.derivative(0)) + std::numeric_limits<real_t>::max() * x * x;
    }
    return c_ * log(chab_ + r);
  }
  real_t c_;
  fvar_t chab_, k_;
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
  boost::math::quadrature::tanh_sinh<real_t> integrator;
  auto logZ = log(real_t(2)) / 2 + 2 * integrator.integrate(func(Jx, Jy, beta), 0, pi / 2);
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
  auto gamma0 = log((1 + cosh(2*a)) / sinh(2*a)) - 2*b;
  for (int_t k = 0; k < 2 * Ly; ++k) {
    auto cosh_g = (cosh(2*a) * cosh(2*b) - cos(pi * k / Ly) * sinh(2*b)) / sinh(2*a);
    auto gamma = (k == 0) ? abs(gamma0) : (log(cosh_g + sqrt(cosh_g * cosh_g - 1)));
    if ((k & 1) == 1) {
      lp0 += (Lx * gamma/2) + log(1 + exp(-(Lx * gamma)));
      lp1 += (Lx * gamma/2) + log(1 - exp(-(Lx * gamma)));
    } else {
      lp2 += (Lx * gamma/2) + log(1 + exp(-(Lx * gamma)));
      lp3 += (Lx * gamma/2) + log(1 - exp(-(Lx * gamma)));
    }
  }
  auto logZ = -log(real_t(2)) / (Lx * Ly) + a + real_t(1)/2 * log(1 - exp(-4*a));
  if (gamma0 > 0) {
    logZ += (lp0 + log(1 + exp(lp1-lp0) + exp(lp2-lp0) - exp(lp3-lp0))) / (Lx * Ly);
  } else if (gamma0 < 0) {
    logZ += (lp0 + log(1 + exp(lp1-lp0) + exp(lp2-lp0) + exp(lp3-lp0))) / (Lx * Ly);
  } else {
    logZ += (lp0 + log(1 + exp(lp1-lp0) + exp(lp2-lp0))) / (Lx * Ly);
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
