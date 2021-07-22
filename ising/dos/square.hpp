/*
   Copyright (C) 2016-2021 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
                              Chihiro Kondo <chihiro.kondo@phys.s.u-tokyo.ac.jp>

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

#pragma once

#include <cmath>
#include <limits>
#include <vector>
#include <boost/math/constants/constants.hpp>
#include <boost/math/differentiation/autodiff.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <lattice/graph.hpp>

namespace ising {
namespace dos {
namespace square {

typedef unsigned long uint_t;

std::vector<uint_t> count(uint_t Lx, uint_t Ly) {
  auto basis = lattice::basis::simple(2);
  auto unitcell = lattice::unitcell(2);
  unitcell.add_site(lattice::coordinate(0, 0), 0);
  unitcell.add_bond(0, 0, lattice::offset(1, 0), 0);
  unitcell.add_bond(0, 0, lattice::offset(0, 1), 0);
  auto graph = lattice::graph(basis, unitcell, lattice::extent(Lx, Ly));
  if (graph.num_bonds() >= std::numeric_limits<uint_t>::digits)
    throw std::range_error("Error: system size is too large\n");
  
  std::vector<uint_t> dos(graph.num_bonds() + 1, 0);
  uint_t num_states = 1 << graph.num_sites();
  for (uint_t c = 0; c < num_states; ++c) {
    uint_t energy = 0;
    for (uint_t b = 0; b < graph.num_bonds(); ++b) {
      uint_t ci = (c >> graph.source(b)) & 1;
      uint_t cj = (c >> graph.target(b)) & 1;
      energy += (ci ^ cj);
    }
    ++dos[energy];
  }
  return dos;
}

// Density of state of square lattice Ising model

// Ref: P. Beale, Phys. Rev. Lett. 76, 78-81 (1996)

namespace {
  
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

}

using namespace boost::math;
using namespace boost::math::differentiation;
namespace mp = boost::multiprecision;
typedef mp::cpp_int int_type;
  
template<unsigned Order, unsigned Digits10, class ExponentType = boost::int32_t>
std::vector<mp::cpp_int> finite(uint_t m, uint_t n) {
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

  std::vector<int_type> dos;
  for (unsigned i = 0; i <= 2 * m * n; ++i) dos.push_back(int_type(zmn.at(i) + 0.01));
  auto sum = std::accumulate(dos.begin(), dos.end(), int_type(0));
  if (sum != power_n(int_type(2), m * n)) {
    std::cerr << "Error: result check failed\n";
    throw(0);
  }
  return dos;
}

}
}
}
