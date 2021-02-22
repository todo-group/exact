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
#include <boost/math/special_functions/factorials.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>

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

template<unsigned Order, unsigned Digits10, class ExponentType = boost::int32_t>
void calc_dos(unsigned m, unsigned n) {
  using namespace boost::math;
  using namespace boost::math::differentiation;
  namespace mp = boost::multiprecision;

  typedef mp::cpp_int int_type;
  typedef mp::number<mp::cpp_dec_float<Digits10, ExponentType>> real_type;
  std::cout << std::setprecision(std::numeric_limits<real_type>::max_digits10);

  auto const x = make_fvar<real_type, Order>(0);
  auto const xpo = make_fvar<real_type, Order>(-1);
  auto const xmo = make_fvar<real_type, Order>(+1);

  auto beta = 2 * x * xpo * xmo;
  auto beta_m = power_n(beta, m);
  auto c0 = power_n(xmo, m) + power_n(x * xpo, m);
  auto s0 = power_n(xmo, m) - power_n(x * xpo, m);
  auto cn = power_n(xpo, m) + power_n(x * xmo, m);
  auto sn = power_n(xpo, m) - power_n(x * xmo, m);

  auto z1 = zero(x);
  auto z2 = zero(x);
  auto z3 = zero(x);
  auto z4 = zero(x);
  if ((n & 1) == 0) {
    z1 = one(x) / 2;
    z2 = one(x) / 2;
    z3 = c0 * cn / 2;
    z4 = s0 * sn / 2;
  } else {
    z1 = cn / 2;
    z2 = sn / 2;
    z3 = c0 / 2;
    z4 = s0 / 2;
  }
  for (unsigned k = 1; k < n; ++k) {
    auto ak = alpha(x, beta, n, k);
    auto v = zero(x);
    for (unsigned j = 0; j <= m; j += 2)
      v += binomial_coefficient<real_type>(m, j)
        * power_n(ak * ak - beta * beta, j/2) * power_n(ak, m-j);
    if ((k & 1) == 1) {
      z1 *= (v + beta_m) / power_n(real_type(2), m - 1);
      z2 *= (v - beta_m) / power_n(real_type(2), m - 1);
    } else {
      z3 *= (v + beta_m) / power_n(real_type(2), m - 1);
      z4 *= (v - beta_m) / power_n(real_type(2), m - 1);
    }
  }
  auto zmn = z1 + z2 + z3 + z4;

  // checking sum rule
  int_type sum;
  for (unsigned i = 0; i <= 2 * m * n; i += 2) sum += int_type(zmn.at(i) + 0.01);
  if (sum != power_n(int_type(2), m * n)) {
    std::cerr << "Error: result check failed\n";
    throw(0);
  }

  for (unsigned i = 0; i <= 2 * m * n; i += 2) std::cout << int_type(zmn.at(i) + 0.01) << ' ';
  std::cout << std::endl;
}

int main(int argc, char **argv) {
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

  if (2 * m * n <= 128) {
    calc_dos<128, 100>(m, n);
  } else if (2 * m * n <= 512) {
    calc_dos<512, 100>(m, n);
  } else if (2 * m * n <= 2048) {
    calc_dos<2048, 100>(m, n);
  } else if (2 * m * n <= 4608) {
    calc_dos<4608, 100>(m, n);
  } else {
    std::cerr << "Error: m * n is too large\n";
    return 127;
  }
}
