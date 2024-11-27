/*
   Copyright (C) 2016-2021 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>

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

#include <iomanip>
#include <gtest/gtest.h>
#include "afh/free_energy/chain_finite.hpp"

TEST(AFHFreeEnergyTest, Chain1) {
  std::cout << std::scientific << std::setprecision(16);
  afh::free_energy::chain::finite solver(6);
  double t = 0.1;
  double eps = 1.0e-10;
  EXPECT_NEAR(-0.467129272955, solver.gs_energy(), eps);
  EXPECT_NEAR( 0.684741648982, solver.gap(), eps);
  double f, e, s;
  std::tie(f, e, s) = solver.free_energy(t);
  EXPECT_NEAR(-0.467182360817, f, eps);
  EXPECT_NEAR(-0.466765889271, e, eps);
  EXPECT_NEAR( 0.004164715458, s, eps);
}
