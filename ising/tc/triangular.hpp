/*****************************************************************************
*
* Copyright (C) 2015-2021 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Critical temperature of triangular lattice Ising model

// reference: J. Stephenson, J. of Math. Phys. 11, 420 (1970)

#pragma once

#include <cmath>
#include <stdexcept>
#include <boost/math/differentiation/autodiff.hpp>
#include <standards/newton.hpp>

namespace {
  
template<typename T>
struct func {
  func(T Ja, T Jb, T Jc) : Ja_(Ja), Jb_(Jb), Jc_(Jc) {}
  auto operator()(T beta) const -> decltype(boost::math::differentiation::make_fvar<T, 1>(beta)) {
    using std::exp;
    auto beta_fvar = boost::math::differentiation::make_fvar<T, 1>(beta);
    auto za = exp(-2 * beta_fvar * Ja_);
    auto zb = exp(-2 * beta_fvar * Jb_);
    auto zc = exp(-2 * beta_fvar * Jc_);
    return za * zb + zb * zc + zc * za - 1;
  }
  T Ja_, Jb_, Jc_;
};

}

namespace ising {
namespace tc {

template<typename T>
inline T triangular(T Ja, T Jb, T Jc) {
  if (Ja * Jb * Jc <= 0) throw(std::invalid_argument("Ja * Jb * Jc should be positive"));
  auto result = standards::newton_1d(func<T>(Ja, Jb, Jc), 1 / (2 * (Ja + Jb + Jc)));
  if (!result.second) throw(std::runtime_error("convergence error"));
  return 1 / result.first;
}

} // end namespace tc
} // end namespace ising
