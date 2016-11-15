/*****************************************************************************
*
* Copyright (C) 1997-2016 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <iostream>
#include "exp_number.hpp"

using namespace log_sum_exp;

int main() {
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  exp_double x = 3;
  exp_double y = 2;
  exp_double z = x + y;
  exp_double v = exp_value(10000);
  exp_double w = exp_value(9999);
  exp_double p = -2.5;

  std::cout << "x = " << x << " = " << static_cast<double>(x) << std::endl;
  std::cout << "y = " << y << " = " << static_cast<double>(y) << std::endl;
  std::cout << "x + y = " << x + y << " = " << static_cast<double>(x + y) << std::endl;
  std::cout << "x - y = " << x - y << " = " << static_cast<double>(x - y) << std::endl;
  std::cout << "x * y = " << x * y << " = " << static_cast<double>(x * y) << std::endl;
  std::cout << "x / y = " << x / y << " = " << static_cast<double>(x / y) << std::endl;
  std::cout << "x + 1.2 = " << x + 1.2 << " = " << static_cast<double>(x + 1.2) << std::endl;
  std::cout << "x - 1.2 = " << x - 1.2 << " = " << static_cast<double>(x - 1.2) << std::endl;
  std::cout << "x - 3.5 = " << x - 3.5 << " = " << static_cast<double>(x - 3.5) << std::endl;
  std::cout << "x * 1.2 = " << x * 1.2 << " = " << static_cast<double>(x * 1.2) << std::endl;
  std::cout << "x * 2 = " << x * 2 << " = " << static_cast<double>(x * 2) << std::endl;
  std::cout << "x / 1.2 = " << x / 1.2 << " = " << static_cast<double>(x / 1.2) << std::endl;
  std::cout << "x / 3 = " << x / 3 << " = " << static_cast<double>(x / 3) << std::endl;
  std::cout << "3.5 + x = " << 3.5 + x << " = " << static_cast<double>(3.5 + x) << std::endl;
  std::cout << "3 + x = " << 3 + x << " = " << static_cast<double>(3 + x) << std::endl;
  std::cout << "3.5 - x = " << 3.5 - x << " = " << static_cast<double>(3.5 - x) << std::endl;
  std::cout << "5 - x = " << 5 - x << " = " << static_cast<double>(4 - x) << std::endl;
  std::cout << "1.2 - x = " << 1.2 - x << " = " << static_cast<double>(1.2 - x) << std::endl;
  std::cout << "3.5 * x = " << 3.5 * x << " = " << static_cast<double>(3.5 * x) << std::endl;
  std::cout << "3 * x = " << 3 * x << " = " << static_cast<double>(3 * x) << std::endl;
  std::cout << "3.5 / x = " << 3.5 / x << " = " << static_cast<double>(3.5 / x) << std::endl;
  std::cout << "3 / x = " << 3 / x << " = " << static_cast<double>(3 / x) << std::endl;

  std::cout << "v = " << v << " = " << static_cast<double>(v) << std::endl;
  std::cout << "w = " << w << " = " << static_cast<double>(w) << std::endl;
  std::cout << "pow(v,3) = " << pow(v,3) << " = " << static_cast<double>(pow(v,3))
            << std::endl;
  std::cout << "v + w = " << v + w << " = " << static_cast<double>(v + w) << std::endl;
  std::cout << "v - w = " << v - w << " = " << static_cast<double>(v - w) << std::endl;
  std::cout << "v * w = " << v * w << " = " << static_cast<double>(v * w) << std::endl;
  std::cout << "v / w = " << v / w << " = " << static_cast<double>(v / w) << std::endl;
  std::cout << "v > w = " << (v > w) << std::endl;
  std::cout << "v >= w = " << (v >= w) << std::endl;
  std::cout << "v == w = " << (v == w) << std::endl;
  std::cout << "v <= w = " << (v <= w) << std::endl;
  std::cout << "v < w = " << (v < w) << std::endl;
  std::cout << "v > v = " << (v > v) << std::endl;
  std::cout << "v >= v = " << (v >= v) << std::endl;
  std::cout << "v == v = " << (v == v) << std::endl;
  std::cout << "v <= v = " << (v <= v) << std::endl;
  std::cout << "v < v = " << (v < v) << std::endl;

  std::cout << "p = " << p << " = " << static_cast<double>(p) << std::endl;
  std::cout << "-p = " << -p << " = " << static_cast<double>(-p) << std::endl;
  std::cout << "p - p = " << p - p << " = " << static_cast<double>(p - p) << std::endl;
  std::cout << "x + p = " << x + p << " = " << static_cast<double>(x + p) << std::endl;
  std::cout << "x - p = " << x - p << " = " << static_cast<double>(x - p) << std::endl;
  std::cout << "x * p = " << x * p << " = " << static_cast<double>(x * p) << std::endl;
  std::cout << "x / p = " << x / p << " = " << static_cast<double>(x / p) << std::endl;
  std::cout << "y + p = " << y + p << " = " << static_cast<double>(y + p) << std::endl;
  std::cout << "y - p = " << y - p << " = " << static_cast<double>(y - p) << std::endl;
  std::cout << "y * p = " << y * p << " = " << static_cast<double>(y * p) << std::endl;
  std::cout << "y / p = " << y / p << " = " << static_cast<double>(y / p) << std::endl;
  std::cout << "x > p = " << (x > p) << std::endl;
  std::cout << "x >= p = " << (x >= p) << std::endl;
  std::cout << "x == p = " << (x == p) << std::endl;
  std::cout << "x <= p = " << (x <= p) << std::endl;
  std::cout << "x < p = " << (x < p) << std::endl;
  std::cout << "y > p = " << (y > p) << std::endl;
  std::cout << "y >= p = " << (y >= p) << std::endl;
  std::cout << "y == p = " << (y == p) << std::endl;
  std::cout << "y <= p = " << (y <= p) << std::endl;
  std::cout << "y < p = " << (y < p) << std::endl;
  std::cout << "p == p = " << (p == p) << std::endl;
  std::cout << "p == -p = " << (p == -p) << std::endl;
  std::cout << "x > 2 = " << (x > 2) << std::endl;
  std::cout << "x >= 2 = " << (x >= 2) << std::endl;
  std::cout << "x == 2 = " << (x == 2) << std::endl;
  std::cout << "x <= 2 = " << (x <= 2) << std::endl;
  std::cout << "x < 2 = " << (x < 2) << std::endl;

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
