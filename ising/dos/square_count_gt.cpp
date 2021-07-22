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

// Density of state of square lattice Ising model

#include <gtest/gtest.h>
#include "square.hpp"

TEST(ising_dos_square, count) {
  auto dos = ising::dos::square::count(4, 4);
  EXPECT_EQ(33, dos.size());
  EXPECT_EQ(2, dos[0]);
  EXPECT_EQ(0, dos[2]);
  EXPECT_EQ(32, dos[4]);
  EXPECT_EQ(64, dos[6]);
  EXPECT_EQ(20524, dos[16]);
  EXPECT_EQ(2, dos[32]);
}
