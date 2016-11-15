// Copyright (C) 2016 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <iostream>
#include "simpson.hpp"

double func1(double x) { return std::sin(x); }
double func2(double x, double y) { return std::sin(x) * std::sin(y); }

int main() {
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  for (unsigned int n = 2; n <= 16; n += 2) {
    std::cout << n << ' '
              << integral::simpson_1d(func1, 0, M_PI, n) - 2 << std::endl;
  }
  if (std::abs(integral::simpson_1d(func1, 0, M_PI, 16) - 2) > 1.0e05) {
    std::cerr << "test failed\n";
    return 1;
  }

  for (unsigned int n = 2; n <= 16; n += 2) {
    std::cout << n << ' '
              << integral::simpson_2d(func2, 0, 0, M_PI, M_PI, n, n) - 4 << std::endl;
  }
  if (std::abs(integral::simpson_2d(func2, 0, 0, M_PI, M_PI, 16, 16) - 4) > 1.0e05) {
    std::cerr << "test failed\n";
    return 1;
  }
  std::cerr << "OK\n";
  
#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& exc) {
  std::cerr << exc.what() << "\n";
  return -1;
}
catch (...) {
  std::cerr << "Fatal Error: Unknown Exception!\n";
  return -2;
}
#endif
  return 0;
}
