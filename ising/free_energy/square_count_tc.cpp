/*****************************************************************************
*
* Copyright (C) 2015-2021 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Free energy, energy, and specific heat of square lattice Ising model

#include <iomanip>
#include <iostream>
#include <boost/math/differentiation/autodiff.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "ising/mp_wrapper.hpp"
#include "square.hpp"

template<typename T, typename I>
void calc(I Lx, I Ly) {
  typedef T real_t;
  namespace ifs = ising::free_energy::square;
  std::cout << std::scientific << std::setprecision(std::numeric_limits<real_t>::digits10);
  real_t beta = log(sqrt(real_t(2)) + 1) / 2;
  auto beta_fvar = boost::math::differentiation::make_fvar<real_t, 2>(beta);
  auto f = ifs::finite_count(Lx, Ly, real_t(1), real_t(1), beta_fvar);
  std::cout << ifs::free_energy(f, beta_fvar) << ' ' << ifs::energy(f, beta_fvar) << ' '
            << ifs::specific_heat(f, beta_fvar) << ' '
            << std::numeric_limits<real_t>::digits10 << std::endl;
}

int main(int argc, char **argv) {
  namespace mp = boost::multiprecision;
  unsigned long Lx, Ly;
  try {
    if (argc == 2) {
      Lx = Ly = std::atoi(argv[1]);
    } else if (argc == 3) {
      Lx = std::atoi(argv[1]);
      Ly = std::atoi(argv[2]);
    } else throw(0);
  } catch (...) {
    std::cerr << "Usage: " << argv[0] << " L\n";
    std::cerr << "       " << argv[0] << " Lx Ly\n";
    return 127;
  }
  calc<float>(Lx, Ly);
  calc<double>(Lx, Ly);
  calc<mp_wrapper<mp::number<mp::cpp_dec_float<30>>>>(Lx, Ly);
  calc<mp_wrapper<mp::cpp_dec_float_50>>(Lx, Ly);
}
