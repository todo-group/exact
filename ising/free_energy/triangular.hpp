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

#pragma once

#include <cmath>
#include <stdexcept>
#include <boost/math/constants/constants.hpp>
#include <boost/math/differentiation/autodiff.hpp>
#include <standards/simpson.hpp>

namespace ising {
namespace free_energy {
namespace triangular {

namespace {

template<typename FVAR, typename T>
struct functor {
  functor(FVAR beta, T Ja, T Jb, T Jc) {
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
functor<FVAR, T> func(FVAR beta, T Ja, T Jb, T Jc) {
  return functor<FVAR, T>(beta, Ja, Jb, Jc);
}

}

template<typename T>
inline std::tuple<T, T, T> infinite(T beta, T Ja, T Jb, T Jc) {
  typedef T real_t;
  using std::log;
  if (beta <= 0) throw(std::invalid_argument("beta should be positive"));
  if (Ja * Jb * Jc <= 0) throw(std::invalid_argument("Ja * Jb * Jc should be positive"));
  real_t pi = boost::math::constants::pi<real_t>();
  auto beta_fvar = boost::math::differentiation::make_fvar<real_t, 2>(beta);
  auto logZ = log(real_t(2)) +
    standards::simpson_2d(func(beta_fvar, Ja, Jb, Jc), real_t(0), real_t(0), 2*pi, 2*pi, 8, 8);
  real_t d2 = logZ.derivative(2);
  for (unsigned long n = 16; n < 1024; n *= 2) {
    logZ = log(real_t(2)) +
      standards::simpson_2d(func(beta_fvar, Ja, Jb, Jc), real_t(0), real_t(0), 2*pi, 2*pi, n, n);
    if (beta * beta * abs((logZ.derivative(2) - d2) / logZ.derivative(0)) <
        10 * std::numeric_limits<real_t>::epsilon()) break;
    d2 = logZ.derivative(2);
  }
  real_t free_energy = - logZ.derivative(0) / beta;
  real_t energy = - logZ.derivative(1);
  real_t specific_heat = beta * beta * logZ.derivative(2);
  return std::make_tuple(free_energy, energy, specific_heat);
}

} // end namespace triangular
} // end namespace free_energy
} // end namespace ising
