/*****************************************************************************
*
* Copyright (C) 2015-2020 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Critical temperature of hexagonal lattice Ising model

#ifndef ISING_TC_HEXAGONAL_HPP
#define ISING_TC_HEXAGONAL_HPP

#include <cmath>
#include <stdexcept>
#include <boost/math/differentiation/autodiff.hpp>
#include <standards/newton.hpp>

namespace {
  
template<typename T>
struct func {
  func(T Ja, T Jb, T Jc) : Ja_(Ja), Jb_(Jb), Jc_(Jc) {}
  auto operator()(T beta) const -> decltype(boost::math::differentiation::make_fvar<T, 2>(beta)) {
    using std::tanh;
    auto beta_fvar = boost::math::differentiation::make_fvar<T, 2>(beta);
    auto za = tanh(beta_fvar * Ja_);
    auto zb = tanh(beta_fvar * Jb_);
    auto zc = tanh(beta_fvar * Jc_);
    return za * zb + zb * zc + zc * za - 1;
  }
  T Ja_, Jb_, Jc_;
};

}

namespace ising {
namespace tc {

template<typename T>
inline T hexagonal(T Ja, T Jb, T Jc) {
  Ja = abs(Ja);
  Jb = abs(Jb);
  Jc = abs(Jc);
  if (Ja * Jb * Jc <= 0) throw(std::invalid_argument("Ja * Jb * Jc should be non-zero"));
  auto result = standards::newton_1d(func<T>(Ja, Jb, Jc), 1 / (Ja + Jb + Jc));
  if (!result.second) throw(std::runtime_error("convergence error"));
  return 1 / result.first;
}

} // end namespace tc
} // end namespace ising

#endif // ISING_TC_HEXAGONAL_HPP
