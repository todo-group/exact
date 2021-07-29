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

// Free energy, energy, and specific heat of triangular lattice Ising model

// reference: R. M. F. Houtappel, Physica 16, 425 (1950)

#pragma once

#include <cmath>
#include <stdexcept>
#include <boost/math/constants/constants.hpp>
#include <boost/math/differentiation/autodiff.hpp>
#include <standards/simpson.hpp>
#include "common.hpp"

namespace ising {
namespace free_energy {
namespace triangular {

namespace {

template<typename T, typename FVAR>
struct functor {
  functor(T Ja, T Jb, T Jc, FVAR beta) {
    using std::cosh; using std::sinh;
    T pi = boost::math::constants::pi<T>();
    c_ = 1 / (8 * pi * pi);
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

template<typename T, typename FVAR>
functor<T, FVAR> func(T Ja, T Jb, T Jc, FVAR beta) {
  return functor<T, FVAR>(Ja, Jb, Jc, beta);
}

}

template<typename T, typename U>
inline U infinite(T Ja, T Jb, T Jc, U beta) {
  const unsigned long max_n = 1 << 16;
  typedef T real_t;
  if (Ja * Jb * Jc <= 0)
    throw(std::invalid_argument("Ja * Jb * Jc should be positive"));
  if (beta <= 0)
    throw(std::invalid_argument("beta should be positive"));
  real_t pi = boost::math::constants::pi<real_t>();
  auto logZ = log(real_t(2)) +
    standards::simpson_2d(func(Ja, Jb, Jc, beta), real_t(0), real_t(0), 2*pi, 2*pi, 8, 8);
  auto pz = logZ;
  auto pe = abs(beta * beta * logZ.derivative(2) / logZ);
  for (unsigned long n = 16; n <= max_n; n *= 2) {
    logZ = log(real_t(2)) +
      standards::simpson_2d(func(Ja, Jb, Jc, beta), real_t(0), real_t(0), 2*pi, 2*pi, n, n);
    auto pn = abs(beta * beta * (logZ - pz).derivative(2) / logZ);
    if (pn < 2 * std::numeric_limits<real_t>::epsilon() || pn > pe) break;
    if (n == max_n) std::cerr << "Warning: integration not converge with n = " << max_n << std::endl;
    pz = logZ;
    pe = pn;
  }
  return - logZ / beta;
}

} // end namespace triangular
} // end namespace free_energy
} // end namespace ising
