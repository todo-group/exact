/*****************************************************************************
*
* Copyright (C) 2015-2021 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <cmath>
#include <random>
#include <gtest/gtest.h>
#include "ising/free_energy/square.hpp"
#include "ising/tc/square.hpp"

namespace ifs = ising::free_energy::square;
namespace autofiff = boost::math::differentiation;
typedef double real_t;

std::mt19937 engine(1234u);
std::uniform_real_distribution<> uniform(0.1, 2.0);

TEST(TcTest, Square1) {
  for (std::size_t i = 0; i < 10; ++i) {
    real_t Jx = uniform(engine);
    real_t Jy = uniform(engine);
    real_t t = uniform(engine) * ising::tc::square(Jx, Jy);
    auto beta = autofiff::make_fvar<real_t, 2>(1 / t);
    auto rf = ising::free_energy::square::finite(4, 2, Jx, Jy, beta);
    auto rc = ising::free_energy::square::finite_count(4, 2, Jx, Jy, beta);
    EXPECT_NEAR(ifs::free_energy(rf, beta), ifs::free_energy(rc, beta), 1e-10);
    EXPECT_NEAR(ifs::energy(rf, beta), ifs::energy(rc, beta), 1e-10);
    EXPECT_NEAR(ifs::specific_heat(rf, beta), ifs::specific_heat(rc, beta), 1e-10);
  }
}

TEST(TcTest, Square2) {
  for (std::size_t i = 0; i < 10; ++i) {
    real_t Jx = uniform(engine);
    real_t Jy = uniform(engine);
    real_t t = uniform(engine) * ising::tc::square(Jx, Jy);
    auto beta = autofiff::make_fvar<real_t, 2>(1 / t);
    auto rf = ising::free_energy::square::finite(4, 4, Jx, Jy, beta);
    auto rc = ising::free_energy::square::finite_count(4, 4, Jx, Jy, beta);
    EXPECT_NEAR(ifs::free_energy(rf, beta), ifs::free_energy(rc, beta), 1e-10);
    EXPECT_NEAR(ifs::energy(rf, beta), ifs::energy(rc, beta), 1e-10);
    EXPECT_NEAR(ifs::specific_heat(rf, beta), ifs::specific_heat(rc, beta), 1e-10);
  }
}

TEST(TcTest, Square3) {
  for (std::size_t i = 0; i < 10; ++i) {
    real_t Jx = uniform(engine);
    real_t Jy = uniform(engine);
    real_t t = uniform(engine) * ising::tc::square(Jx, Jy);
    auto beta = autofiff::make_fvar<real_t, 2>(1 / t);
    auto rf = ising::free_energy::square::finite(10000, 10000, Jx, Jy, beta);
    auto ri = ising::free_energy::square::infinite(Jx, Jy, beta);
    EXPECT_NEAR(ifs::free_energy(rf, beta), ifs::free_energy(ri, beta), 1e-6);
    EXPECT_NEAR(ifs::energy(rf, beta), ifs::energy(ri, beta), 1e-10);
    EXPECT_NEAR(ifs::specific_heat(rf, beta), ifs::specific_heat(ri, beta), 1e-10);
  }
}
