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

// Critical temperature of square lattice Ising model

// reference: L. Onsager, Phys. Rev. 65, 117 (1944)

#pragma once

#include <cmath>
#include <stdexcept>
#include <boost/math/differentiation/autodiff.hpp>
#include <standards/newton.hpp>

namespace {
  
template<typename T>
struct func {
  func(T Jx, T Jy) : Jx_(Jx), Jy_(Jy) {}
  auto operator()(T beta) const -> decltype(boost::math::differentiation::make_fvar<T, 1>(beta)) {
    using std::sinh;
    auto beta_fvar = boost::math::differentiation::make_fvar<T, 1>(beta);
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
