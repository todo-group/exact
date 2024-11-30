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
#include "ising/free_energy/square_tm.hpp"
#include "common.hpp"

using namespace ising::free_energy;

TEST(IsingFreeEnergy, SquareCount0) {
  unsigned Lx = 4;
  unsigned Ly = 4;
  double Jx = 1.5;
  double Jy = 2.5;
  double t = 2;
  double h = 0;
  square::transfer_matrix<double>::result_t beta;
  beta.set(0, 0, 1 / t);
  auto f = square::transfer_matrix<double>::calc(Lx, Ly, Jx, Jy, 1 / t, h);
  EXPECT_NEAR(-4.087359662653047e+00, free_energy(f, beta, h), 1e-10); 
  EXPECT_NEAR(-3.994108759068211e+00, energy(f, beta, h), 1e-10);
  EXPECT_NEAR(2.452622208849045e-02, specific_heat(f, beta, h), 1e-10);
  EXPECT_NEAR(1.597700713244840e+01, magnetization2(f, beta, h), 1e-10);
}
