/*****************************************************************************
*
* Copyright (C) 2015-2021 by Synge Todo <wistaria@phy.s.u-tokyo.ac.jp>
*                            Chihiro Kondo <chihiro.kondo@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Density of state of square lattice Ising model

// Ref: P. Beale, Phys. Rev. Lett. 76, 78-81 (1996)

#include <cmath>
#include <boost/math/constants/constants.hpp>
#include <boost/math/differentiation/autodiff.hpp>
#include <boost/math/special_functions/binomial.hpp>

template<typename T>
T zero(const T& x) { return 0 * x; }

template<typename T>
T one(const T& x) { return 1 + zero(x); }

template<typename T>
T power_n(const T& x, unsigned n) {
  auto res = one(x);
  for (unsigned i = 0; i < n; ++i) res *= x;
  return res;
}

template<typename T>
T alpha(const T& x, const T& beta, unsigned n, unsigned k) {
  using std::cos;
  auto pi = boost::math::constants::pi<typename T::root_type>();
  return power_n(1 + x * x, 2) - beta * cos(pi * k / n);
}

template<typename T>
T c2(const T& x, const T& beta, unsigned m, unsigned n, unsigned k) {
  using namespace boost::math;
  auto ak = alpha(x, beta, n, k);
  auto res = zero(x);
  for (unsigned j = 0; j <= m; j += 2)
    res += binomial_coefficient<typename T::root_type>(m, j)
      * power_n(ak * ak - beta * beta, j/2) * power_n(ak, m-j);
  res += power_n(beta, m);
  return res / power_n(2, m - 1);
}

template<typename T>
T s2(const T& x, const T& beta, unsigned m, unsigned n, unsigned k) {
  using namespace boost::math;
  auto ak = alpha(x, beta, n, k);
  auto res = zero(x);
  for (unsigned j = 0; j <= m; j += 2)
    res += binomial_coefficient<typename T::root_type>(m, j)
      * power_n(ak * ak - beta * beta, j/2) * power_n(ak, m-j);
  res -= power_n(beta, m);
  return res / power_n(2, m - 1);
}

int main(int argc, char **argv) {
  using namespace boost::math;
  using namespace boost::math::differentiation;

  constexpr unsigned ORDER = 128;
  unsigned m, n;
  try {
    if (argc == 2) {
      m = n = std::stoi(argv[1]);
    } else if (argc == 3) {
      m = std::stoi(argv[1]);
      n = std::stoi(argv[2]);
    } else throw(0);
    if (m == 0 || n == 0) {
      std::cerr << "Error: Lx and Ly should be positive\n";
      throw(0);
    }
    if ((m & 1) == 1) {
      std::cerr << "Error: Lx should be even\n";
      throw(0);
    }
  } catch (...) {
    std::cerr << "Usage: " << argv[0] << " L\n";
    std::cerr << "       " << argv[0] << " Lx Ly\n";
    return 127;
  }

  auto const x = make_fvar<double, ORDER>(0.0);

  auto beta = 2 * x * (1 - x * x);
  auto c0 = power_n(1 - x, m) + power_n(x * (1 + x), m);
  auto s0 = power_n(1 - x, m) - power_n(x * (1 + x), m);
  auto cn = power_n(1 + x, m) + power_n(x * (1 - x), m);
  auto sn = power_n(1 + x, m) - power_n(x * (1 - x), m);

  auto z1 = zero(x);
  auto z2 = zero(x);
  auto z3 = zero(x);
  auto z4 = zero(x);
  if ((n & 1) == 0) {
    z1 = z2 = one(x) / 2;
    z3 = c0 * cn / 2;
    z4 = s0 * sn / 2;
  } else {
    z1 = cn / 2;
    z2 = sn / 2;
    z3 = c0 / 2;
    z4 = s0 / 2;
  }
  for (unsigned k = 1; k < n; ++k) {
    if ((k & 1) == 1) {
      z1 *= c2(x, beta, m, n, k);
      z2 *= s2(x, beta, m, n, k);
    } else {
      z3 *= c2(x, beta, m, n, k);
      z4 *= s2(x, beta, m, n, k);
    }
  }
  auto zmn = z1 + z2 + z3 + z4;
  for (unsigned i = 0; i <= 2 * m * n; i += 2) std::cout << unsigned(zmn.at(i) + 0.5) << ' ';
  std::cout << std::endl;
}
