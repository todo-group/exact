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
#include <boost/math/constants/constants.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>

namespace tfi {
namespace energy {
namespace chain {

namespace {

template<typename T>
struct functor {
  typedef T real_t;
  functor(real_t J, real_t Gamma) : J_(J), Gamma_(Gamma) {}
  real_t operator()(real_t k) const {
    using std::abs; using std::cos; using std::sin; using std::sqrt; using std::pow;
    return 2 * abs(J_) * sqrt(pow(cos(k) - (Gamma_ / J_), 2) + pow(sin(k), 2));
  }
  real_t J_, Gamma_;
};

template<typename T>
functor<T> func(T J, T Gamma) { return functor<T>(J, Gamma); }

}
  
template<typename T>
inline T infinite(T J, T Gamma) {
  using std::abs;
  typedef T real_t;
  if (abs(J)> 0) {
    auto pi = boost::math::constants::pi<real_t>();
    boost::math::quadrature::tanh_sinh<real_t> integrator;
    return -integrator.integrate(func(J, Gamma), 0, pi) / (2 * pi);
  } else {
    return -abs(Gamma);
  }
}

} // end namespace chain
} // end namespace energy
} // end namespace tfi
