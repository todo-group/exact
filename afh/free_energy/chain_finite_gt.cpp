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

#include <gtest/gtest.h>
#include "afh/free_energy/chain_finite.hpp"

TEST(AFHFreeEnergyTest, Chain1) {
  afh::free_energy::chain::finite solver(6);
  double t = 0.1;
  EXPECT_DOUBLE_EQ(-0.46712927295533307, solver.gs_energy());
  EXPECT_DOUBLE_EQ(0.68474164898209855, solver.gap());
  double f, e, s;
  std::tie(f, e, s) = solver.free_energy(t);
  EXPECT_DOUBLE_EQ(-0.37133552900632999, f);
  EXPECT_DOUBLE_EQ(-0.35286500691579348, e);
  EXPECT_DOUBLE_EQ(0.18470522090536512, s);
}
