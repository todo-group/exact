/*
   Copyright (C) 2015-2020 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>

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
#include <random>
#include <gtest/gtest.h>
#include "ising/free_energy/square.hpp"
#include "ising/tc/square.hpp"

std::mt19937 engine(1234u);
std::uniform_real_distribution<> uniform(0.1, 2.0);

TEST(TcTest, Square1) {
  for (std::size_t i = 0; i < 10; ++i) {
    double Jx = uniform(engine);
    double Jy = uniform(engine);
    double T = uniform(engine) * ising::tc::square(Jx, Jy);
    auto rf = ising::free_energy::square::finite(4, 2, Jx, Jy, T);
    auto rc = ising::free_energy::square::finite_count(4, 2, Jx, Jy, T);
    EXPECT_NEAR(std::get<0>(rf), std::get<0>(rc), 1e-10);
    EXPECT_NEAR(std::get<1>(rf), std::get<1>(rc), 1e-10);
    EXPECT_NEAR(std::get<2>(rf), std::get<2>(rc), 1e-10);
  }
}

TEST(TcTest, Square2) {
  for (std::size_t i = 0; i < 10; ++i) {
    double Jx = uniform(engine);
    double Jy = uniform(engine);
    double T = uniform(engine) * ising::tc::square(Jx, Jy);
    auto rf = ising::free_energy::square::finite(4, 4, Jx, Jy, T);
    auto rc = ising::free_energy::square::finite_count(4, 4, Jx, Jy, T);
    EXPECT_NEAR(std::get<0>(rf), std::get<0>(rc), 1e-10);
    EXPECT_NEAR(std::get<1>(rf), std::get<1>(rc), 1e-10);
    EXPECT_NEAR(std::get<2>(rf), std::get<2>(rc), 1e-10);
  }
}

TEST(TcTest, Square3) {
  for (std::size_t i = 0; i < 10; ++i) {
    double Jx = uniform(engine);
    double Jy = uniform(engine);
    double T = uniform(engine) * ising::tc::square(Jx, Jy);
    auto rf = ising::free_energy::square::infinite(Jx, Jy, T);
    auto rc = ising::free_energy::square::finite(10000, 10000, Jx, Jy, T);
    EXPECT_NEAR(std::get<0>(rf), std::get<0>(rc), 1e-7);
    EXPECT_NEAR(std::get<1>(rf), std::get<1>(rc), 1e-6);
    EXPECT_NEAR(std::get<2>(rf), std::get<2>(rc), 1e-4);
  }
}
