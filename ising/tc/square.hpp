/*****************************************************************************
*
* Copyright (C) 2015-2020 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Critical temperature of square lattice Ising model

// reference: L. Onsager, Phys. Rev. 65, 117 (1944)

#ifndef ISING_TC_SQUARE_HPP
#define ISING_TC_SQUARE_HPP

#include <cmath>
#include <stdexcept>
#include <boost/math/differentiation/autodiff.hpp>
#include <standards/newton.hpp>

namespace {
  
template<typename T>
struct func {
  func(T Jx, T Jy) : Jx_(Jx), Jy_(Jy) {}
  auto operator()(T beta) const -> decltype(boost::math::differentiation::make_fvar<T, 2>(beta)) {
    using std::sinh;
    auto beta_fvar = boost::math::differentiation::make_fvar<T, 2>(beta);
    return sinh(2 * Jx_ * beta_fvar) * sinh(2 * Jy_ * beta_fvar) - 1;
  }
  T Jx_, Jy_;
};

}

namespace ising {
namespace tc {

template<typename T>
inline T square(T Jx, T Jy) {
  Jx = abs(Jx);
  Jy = abs(Jy);
  if (Jx * Jy <= 0) throw(std::invalid_argument("Jx * Jy should be non-zero"));
  auto result = standards::newton_1d(func<T>(Jx, Jy), 1 / (2 * (Jx + Jy)));
  if (!result.second) throw(std::runtime_error("convergence error"));
  return 1 / result.first;
}

} // end namespace tc
} // end namespace ising

#endif // ISING_TC_SQUARE_HPP
