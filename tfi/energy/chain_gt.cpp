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
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "chain.hpp"
#include "convert.hpp"

using namespace boost::multiprecision;
using namespace tfi::energy;

TEST(TFIEnergy, ChainF0) {
  typedef float real_t;
  auto J = convert<real_t>("1");
  auto Gamma = convert<real_t>("1");
  auto ene = chain::infinite(J, Gamma);
  EXPECT_FLOAT_EQ(-1.2732395447351626861510701069801148962756771659236515899813387524711743810738122807209104221300246876, ene);
}

TEST(TFIEnergy, ChainF1) {
  typedef float real_t;
  auto J = convert<real_t>("1.2");
  auto Gamma = convert<real_t>("-0.8");
  auto ene = chain::infinite(J, Gamma);
  EXPECT_FLOAT_EQ(-1.3375409772289557122297861376619311247381065311346132244923834297492659680725295510811206803209285916, ene);
}

TEST(TFIEnergy, ChainD) {
  typedef double real_t;
  auto J = convert<real_t>("1");
  auto Gamma = convert<real_t>("1");
  auto ene = chain::infinite(J, Gamma);
  EXPECT_DOUBLE_EQ(-1.2732395447351626861510701069801148962756771659236515899813387524711743810738122807209104221300246876, ene);
}

TEST(TFIEnergy, ChainD1) {
  typedef double real_t;
  auto J = convert<real_t>("1.2");
  auto Gamma = convert<real_t>("-0.8");
  auto ene = chain::infinite(J, Gamma);
  EXPECT_FLOAT_EQ(-1.3375409772289557122297861376619311247381065311346132244923834297492659680725295510811206803209285916, ene);
}
