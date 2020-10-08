// Copyright (C) 2015-2020 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef STANDARDS_NEWTON_HPP
#define STANDARDS_NEWTON_HPP

#include <cmath>
#include <limits>

namespace standards {

// Solve f(x) = 0
//   Return: pair of solution and number of iteractions
//           iteration = 0 means that algorithm does not converge during max_iter.

template<typename F, typename T>
std::pair<T, std::size_t> newton_1d(F const& func, T x, std::size_t max_iter = 256) {
  using std::abs;
  std::size_t iter = 0;
  auto v = func(x);
  for (; iter < max_iter; ++iter) {
    auto xn = x - v.derivative(0) / v.derivative(1);
    auto vn = func(xn);
    if (abs(xn - x) < 2 * abs(x) * std::numeric_limits<T>::epsilon() ||
        abs(vn.derivative(0)) < 2 * std::numeric_limits<T>::epsilon())
      return std::make_pair(xn, iter + 1);
    x = xn;
    v = vn;
  }
  return std::make_pair(x, 0);
}

}

#endif // STANDARDS_NEWTON_HPP
