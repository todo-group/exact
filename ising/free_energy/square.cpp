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
#include "square.hpp"

int main(int argc, char **argv) {
  namespace ifs = ising::free_energy::square;
  namespace autofiff = boost::math::differentiation;
  typedef double real_t;
  real_t Jx, Jy, tmin, tmax;
  real_t dt = 1;
  try {
    if (argc == 2) {
      Jx = Jy = 1;
      tmin = tmax = std::atof(argv[1]);
    } else if (argc == 3) {
      Jx = Jy = std::atof(argv[1]);
      tmin = tmax = std::atof(argv[2]);
    } else if (argc == 4) {
      Jx = std::atof(argv[1]);
      Jy = std::atof(argv[2]);
      tmin = tmax = std::atof(argv[3]);
    } else if (argc == 6) {
      Jx = std::atof(argv[1]);
      Jy = std::atof(argv[2]);
      tmin = std::atof(argv[3]);
      tmax = std::atof(argv[4]);
      dt = std::atof(argv[5]);
    } else throw(0);
  } catch (...) {
    std::cerr << "Usage: " << argv[0] << " T\n";
    std::cerr << "       " << argv[0] << " J T\n";
    std::cerr << "       " << argv[0] << " Jx Jy T\n";
    std::cerr << "       " << argv[0] << " Jx Jy Tmin Tmax dT\n";
    return 127;
  }
  std::cout << std::scientific << std::setprecision(std::numeric_limits<real_t>::digits10);
  for (real_t t = tmin; t < tmax + 0.01 * dt; t += dt) {
    auto beta = autofiff::make_fvar<real_t, 2>(1 / t);
    auto f = ifs::infinite(Jx, Jy, beta);
    std::cout << Jx << ' ' << Jy << ' ' << t << ' '
              << ifs::free_energy(f, beta) << ' ' << ifs::energy(f, beta) << ' '
              << ifs::specific_heat(f, beta) << std::endl;
  }
}
