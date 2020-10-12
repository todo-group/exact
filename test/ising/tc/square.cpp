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
#include <gtest/gtest.h>
#include "ising/tc/square.hpp"

TEST(TcTest, Square1) {
  double tc = 2 / log(1 + sqrt(2.0));
  EXPECT_DOUBLE_EQ(tc, ising::tc::square(1.0, 1.0));
}

