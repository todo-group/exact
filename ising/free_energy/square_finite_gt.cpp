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

TEST(IsingFreeEnergy, SquareFinite0) {
  typedef double real_t;
  unsigned Lx = 4;
  unsigned Ly = 4;
  auto Jx = convert<real_t>("1.5");
  auto Jy = convert<real_t>("2.5");
  auto t = convert<real_t>("2");
  auto vars = boost::math::differentiation::make_ftuple<real_t, 2, 2>(1 / t, 0);
  auto& beta = std::get<0>(vars);
  auto& h = std::get<1>(vars);
  auto fc = square::finite_count(Lx, Ly, Jx, Jy, beta, h);
  auto ff = square::finite(Lx, Ly, Jx, Jy, beta);
  double eps = 1e-12;
  EXPECT_TRUE(abs(free_energy(fc, beta, h) - free_energy(ff, beta)) < eps);
  EXPECT_TRUE(abs(energy(fc, beta, h) - energy(ff, beta)) < eps);
  EXPECT_TRUE(abs(specific_heat(fc, beta, h) - specific_heat(ff, beta)) < eps);
}

TEST(IsingFreeEnergy, SquareFinite1) {
  typedef mp_wrapper<cpp_dec_float_50> real_t;
  unsigned Lx = 4;
  unsigned Ly = 4;
  auto Jx = convert<real_t>("1.5");
  auto Jy = convert<real_t>("2.5");
  auto t = convert<real_t>("2");
  auto vars = boost::math::differentiation::make_ftuple<real_t, 2, 2>(1 / t, 0);
  auto& beta = std::get<0>(vars);
  auto& h = std::get<1>(vars);
  auto fc = square::finite_count(Lx, Ly, Jx, Jy, beta, h);
  auto ff = square::finite(Lx, Ly, Jx, Jy, beta);
  double eps = 1e-40;
  EXPECT_TRUE(abs(free_energy(fc, beta, h) - free_energy(ff, beta)) < eps);
  EXPECT_TRUE(abs(energy(fc, beta, h) - energy(ff, beta)) < eps);
  EXPECT_TRUE(abs(specific_heat(fc, beta, h) - specific_heat(ff, beta)) < eps);
}

TEST(IsingFreeEnergy, SquareFinite2) {
  typedef mp_wrapper<cpp_dec_float_100> real_t;
  unsigned Lx = 4;
  unsigned Ly = 4;
  auto Jx = convert<real_t>("1.5");
  auto Jy = convert<real_t>("2.5");
  auto t = convert<real_t>("2");
  auto vars = boost::math::differentiation::make_ftuple<real_t, 2, 2>(1 / t, 0);
  auto& beta = std::get<0>(vars);
  auto& h = std::get<1>(vars);
  auto fc = square::finite_count(Lx, Ly, Jx, Jy, beta, h);
  auto ff = square::finite(Lx, Ly, Jx, Jy, beta);
  double eps = 1e-80;
  EXPECT_TRUE(abs(free_energy(fc, beta, h) - free_energy(ff, beta)) < eps);
  EXPECT_TRUE(abs(energy(fc, beta, h) - energy(ff, beta)) < eps);
  EXPECT_TRUE(abs(specific_heat(fc, beta, h) - specific_heat(ff, beta)) < eps);
}
