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
  using std::log;
  if (Jx <= 0 || Jy <= 0) throw(std::invalid_argument("Jx and Jy should be positive"));
  if (beta <= 0) throw(std::invalid_argument("beta should be positive"));
  const real_t pi = boost::math::constants::pi<real_t>();
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
  typedef T real_t;
  using std::exp; using std::log;
  if (Lx <= 0 || Ly <= 0) throw(std::invalid_argument("Lx and Ly should be positive"));
  if (Jx <= 0 || Jy <= 0) throw(std::invalid_argument("Jx and Jy should be positive"));
  if (beta <= 0) throw(std::invalid_argument("beta should be positive"));
  const real_t pi = boost::math::constants::pi<real_t>();
  auto beta_fvar = boost::math::differentiation::make_fvar<real_t, 2>(beta);
  auto a = beta_fvar * Jx;
  auto b = beta_fvar * Jy;
  auto a_bar = log((1 + cosh(2*a)) / sinh(2*a)) / 2;
  auto gamma_max = 2 * (a_bar + b);
  decltype(beta_fvar) p0(1), p1(1), p2(1), p3(1);
  for (std::size_t k = 0; k < 2 * Lx; ++k) {
    auto cosh_g = (cosh(2*a) * cosh(2*b) - cos(pi*k/Lx) * sinh(2*b)) / sinh(2*a);
    auto gamma = log(cosh_g + sqrt(cosh_g * cosh_g - 1));
    if (k == 0) gamma = 2 * (a_bar - b);
    if ((k & 1) == 0) {
      p2 *= exp(-(Ly * (gamma_max - gamma) / 2)) + exp(-(Ly * (gamma_max + gamma) / 2));
      p3 *= exp(-(Ly * (gamma_max - gamma) / 2)) - exp(-(Ly * (gamma_max + gamma) / 2));
    } else {
      p0 *= exp(-(Ly * (gamma_max - gamma) / 2)) + exp(-(Ly * (gamma_max + gamma) / 2));
      p1 *= exp(-(Ly * (gamma_max - gamma) / 2)) - exp(-(Ly * (gamma_max + gamma) / 2));
    }
  }
  std::cerr << a.derivative(0) << ' ';
  std::cerr << b.derivative(0) << ' ';
  std::cerr << a_bar.derivative(0) << ' ';
  std::cerr << gamma_max.derivative(0) << ' ';
  std::cerr << p0.derivative(0) << ' ';
  std::cerr << p1.derivative(0) << ' ';
  std::cerr << p2.derivative(0) << ' ';
  std::cerr << p3.derivative(0) << ' ';
  std::cerr << std::endl;
  auto logZ = - log(real_t(2)) / (Lx * Ly) + a + 0.5 * log(1 - exp(-4*a)) + 0.5 * gamma_max + log(p0 + p1 + p2 - p3) / (Lx * Ly);
  real_t free_energy = - logZ.derivative(0) / beta;
  real_t energy = - logZ.derivative(1);
  real_t specific_heat = beta * beta * logZ.derivative(2);
  return std::make_tuple(free_energy, energy, specific_heat);
}

} // end namespace square
} // end namespace free_energy
} // end namespace ising
