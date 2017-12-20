// Copyright (C) 2015-2016 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef INTEGRAL_SIMPSON_HPP
#define INTEGRAL_SIMPSON_HPP

#include <boost/throw_exception.hpp>
#include <cmath>
#include <stdexcept>

namespace integral {

template<typename F>
double simpson_1d(F const& func, double x0, double x1, unsigned int n) {
  if (n == 0 || (n % 2) != 0)
    boost::throw_exception(std::invalid_argument("n should be positive and a multiple of two"));
  double dx = (x1 - x0) / n;
  double g = 0.0;
  // i == 0 || i == n
  g += func(x0) + func(x1);
  // i = 1/2 ... n-1/2
  for (unsigned int i = 0; i < n; ++i) {
    double x = x0 + dx * (i + 0.5);
    g += 4 * func(x);
  }
  // i = 1 ... n-1
  for (unsigned int i = 1; i < n; ++i) {
    double x = x0 + dx * i;
    g += 2 * func(x);
  }
  g *= dx / 6;
  return g;
}

template<typename F>
double simpson_2d(F const& func, double x0, double y0, double x1, double y1,
  unsigned int nx, unsigned int ny) {
  if (nx == 0 || (nx % 2) != 0 || ny == 0 || (ny % 2) != 0)
    boost::throw_exception(std::invalid_argument("nx and ny should be positive and a multiple of two"));
  double dx = (x1 - x0) / nx;
  double dy = (y1 - y0) / ny;
  double g = 0.0;
  // i == 0
  {
    double x = x0;
    //   j == 0 || j == ny
    g += func(x, y0) + func(x, y1);
    //   j = 1/2 ... ny-1/2
    for (unsigned int j = 0; j < ny; ++j) {
      double y = y0 + dy * (j + 0.5);
      g += 4 * func(x, y);
    }
    //   j = 1 ... ny-1
    for (unsigned int j = 1; j < ny; ++j) {
      double y = y0 + dy * j;
      g += 2 * func(x, y);
    }
  }
  // i = 1/2 ... nx-1/2
  for (unsigned int i = 0; i < nx; ++i) {
    double x = x0 + dx * (i + 0.5);
    //   j == 0 || j == ny
    g += 4 * (func(x, y0) + func(x, y1));
    //   j = 1/2 ... ny-1/2
    for (unsigned int j = 0; j < ny; ++j) {
      double y = y0 + dy * (j + 0.5);
      g += 16 * func(x, y);
    }
    //   j = 1 ... ny-1
    for (unsigned int j = 1; j < ny; ++j) {
      double y = y0 + dy * j;
      g += 8 * func(x, y);
    }
  }
  // i = 1 ... nx-1
  for (unsigned int i = 1; i < nx; ++i) {
    double x = x0 + dx * i;
    //   j == 0 || j == ny
    g += 2 * (func(x, y0) + func(x, y1));
    //   j = 1/2 ... ny-1/2
    for (unsigned int j = 0; j < ny; ++j) {
      double y = y0 + dy * (j + 0.5);
      g += 8 * func(x, y);
    }
    //   j = 1 ... ny-1
    for (unsigned int j = 1; j < ny; ++j) {
      double y = y0 + dy * j;
      g += 4 * func(x, y);
    }
  }
  // i == nx
  {
    double x = x1;
    //   j == 0 || j == ny
    g += func(x, y0) + func(x, y1);
    //   j = 1/2 ... ny-1/2
    for (unsigned int j = 0; j < ny; ++j) {
      double y = y0 + dy * (j + 0.5);
      g += 4 * func(x, y);
    }
    //   j = 1 ... ny-1
    for (unsigned int j = 1; j < ny; ++j) {
      double y = y0 + dy * j;
      g += 2 * func(x, y);
    }
  }
  g *= dx * dy / 36;
  return g;
}

} // end namespace integral

#endif // INTEGRAL_SIMPSON_HPP
