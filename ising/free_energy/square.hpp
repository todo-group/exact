/*****************************************************************************
*
* Copyright (C) 2015-2020 by Synge Todo <wistaria@phy.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Calculating free energy, energy, and specific heat of square lattice Ising model

// reference: https://en.wikipedia.org/wiki/Square-lattice_Ising_model

#ifndef ISING_SQUARE_INFINITE_HPP
#define ISING_SQUARE_INFINITE_HPP

#include <cmath>
#include <stdexcept>
#include <boost/math/constants/constants.hpp>
#include <boost/math/differentiation/autodiff.hpp>
#include <standards/simpson.hpp>

namespace {

template<typename FVAR, typename T>
struct func_t {
  func_t(FVAR beta, T Jx, T Jy) {
    using std::cosh; using std::sinh;
    c_ = 1 / (2 * boost::math::constants::pi<T>());
    chab_ = cosh(2 * beta * Jx) * cosh(2 * beta * Jy);
    k_ = 1.0 / (sinh(2 * beta * Jx) * sinh(2 * beta * Jy));
  }
  FVAR operator()(T t) const {
    using std::cos; using std::log;
    return c_ * log(chab_ + sqrt(1 + k_ * k_ - 2 * k_ * cos(2 * t)) / k_);
  }
  T c_;
  FVAR chab_, k_;
};

template<typename FVAR, typename T>
func_t<FVAR, T> func(FVAR beta, T Jx, T Jy) {
  return func_t<FVAR, T>(beta, Jx, Jy);
}
  
}

namespace ising {
namespace square {

template<typename T>
inline std::tuple<T, T, T> infinite(T beta, T Jx, T Jy) {
  using std::log;
  if (beta <= 0) throw(std::invalid_argument("beta should be positive"));
  if (Jx * Jy <= 0) throw(std::invalid_argument("Jx * Jy should be positive"));
  T pi = boost::math::constants::pi<T>();
  auto beta_fvar = boost::math::differentiation::make_fvar<T, 2>(beta);
  auto logZ = log(T(2)) / 2 + standards::simpson_1d(func(beta_fvar, Jx, Jy), T(0), pi, 8);
  T d2 = logZ.derivative(2);
  for (unsigned long n = 16; n < 1024; n *= 2) {
    logZ = log(T(2)) / 2 + standards::simpson_1d(func(beta_fvar, Jx, Jy), T(0), pi, n);
    if (beta * beta * abs((logZ.derivative(2) - d2) / logZ.derivative(0)) <
        10 * std::numeric_limits<T>::epsilon()) break;
    d2 = logZ.derivative(2);
  }
  T free_energy = - logZ.derivative(0) / beta;
  T energy = - logZ.derivative(1);
  T specific_heat = beta * beta * logZ.derivative(2);
  return std::make_tuple(free_energy, energy, specific_heat);
}

} // end namespace square
} // end namespace ising

#endif // ISING_SQUARE_INFINITE_HPP
