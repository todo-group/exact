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

#include <cmath>
#include <gtest/gtest.h>
#include <boost/math/differentiation/autodiff.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "ising/mp_wrapper.hpp"
#include "ising/free_energy/common.hpp"
#include "ising/free_energy/square.hpp"

using namespace boost::multiprecision;
using namespace ising::free_energy;

TEST(IsingFreeEnergy, SquareCount0) {
  typedef double real_t;
  unsigned Lx = 4;
  unsigned Ly = 4;
  auto Jx = convert<real_t>("1.5");
  auto Jy = convert<real_t>("2.5");
  auto t = convert<real_t>("2");
  auto beta = boost::math::differentiation::make_fvar<real_t, 2>(1 / t);
  auto f = square::finite_count(Lx, Ly, Jx, Jy, beta);
  EXPECT_DOUBLE_EQ(-4.087359662653047e+00, free_energy(f, beta));
  EXPECT_DOUBLE_EQ(-3.994108759068211e+00, energy(f, beta));
  EXPECT_DOUBLE_EQ(2.452622208849045e-02, specific_heat(f, beta));
}
