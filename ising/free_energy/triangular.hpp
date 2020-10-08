/*****************************************************************************
*
* Copyright (C) 2015-2020 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Calculating free energy, energy, and specific heat of triangular lattice Ising model

// reference: R. M. F. Houtappel, Physica 16, 425 (1950)

#ifndef ISING_TRIANGULAR_INFINITE_HPP
#define ISING_TRIANGULAR_INFINITE_HPP

#include <cmath>
#include <stdexcept>
#include <boost/math/constants/constants.hpp>
#include <boost/math/differentiation/autodiff.hpp>
#include <standards/simpson.hpp>

namespace {

template<typename FVAR, typename T>
struct func_t {
  func_t(FVAR beta, T Ja, T Jb, T Jc) {
    using std::cosh; using std::sinh;
    T pi = boost::math::constants::pi<T>();
    c_ = 1.0 / (8 * pi * pi);
    sh2a_ = sinh(2 * beta * Ja);
    sh2b_ = sinh(2 * beta * Jb);
    sh2c_ = sinh(2 * beta * Jc);
    cshabc_ = cosh(2 * beta * Ja) * cosh(2 * beta * Jb) * cosh(2 * beta * Jc)
      + sh2a_ * sh2b_ * sh2c_;
  }
  FVAR operator()(T t1, T t2) const {
    using std::cos; using std::log;
    return c_ * log(cshabc_ - sh2a_ * cos(t1) - sh2b_ * cos(t2) - sh2c_ * cos(t1 + t2));
  }
  T c_;
  FVAR sh2a_, sh2b_, sh2c_, cshabc_;
};

template<typename FVAR, typename T>
func_t<FVAR, T> func(FVAR beta, T Ja, T Jb, T Jc) {
  return func_t<FVAR, T>(beta, Ja, Jb, Jc);
}

}

namespace ising {
namespace triangular {

template<typename T>
inline std::tuple<T, T, T> infinite(T beta, T Ja, T Jb, T Jc) {
  using std::log;
  if (beta <= 0) throw(std::invalid_argument("beta should be positive"));
  if (Ja * Jb * Jc <= 0) throw(std::invalid_argument("Ja * Jb * Jc should be positive"));
  T pi = boost::math::constants::pi<T>();
  auto beta_fvar = boost::math::differentiation::make_fvar<T, 2>(beta);
  auto logZ = log(T(2)) +
    standards::simpson_2d(func(beta_fvar, Ja, Jb, Jc), T(0), T(0), 2*pi, 2*pi, 8, 8);
  T d2 = logZ.derivative(2);
  for (unsigned long n = 16; n < 1024; n *= 2) {
    logZ = log(T(2)) +
      standards::simpson_2d(func(beta_fvar, Ja, Jb, Jc), T(0), T(0), 2*pi, 2*pi, n, n);
    if (beta * beta * abs((logZ.derivative(2) - d2) / logZ.derivative(0)) <
        10 * std::numeric_limits<T>::epsilon()) break;
    d2 = logZ.derivative(2);
  }
  T free_energy = - logZ.derivative(0) / beta;
  T energy = - logZ.derivative(1);
  T specific_heat = beta * beta * logZ.derivative(2);
  return std::make_tuple(free_energy, energy, specific_heat);
}

} // end namespace triangular
} // end namespace ising

#endif // ISING_TRIANGULAR_INFINITE_HPP
