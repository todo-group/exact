/*****************************************************************************
*
* Copyright (C) 2015-2020 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Free energy, energy, and specific heat of square lattice Ising model

#include <iomanip>
#include <iostream>
#include <boost/lexical_cast.hpp>
#include "square.hpp"

int main(int argc, char **argv) {
  typedef unsigned long uint_t;
  typedef double real_t;
  uint_t Lx, Ly;
  real_t Jx, Jy, tmin, tmax;
  real_t dt = 1;
  try {
    if (argc == 3) {
      Jx = Jy = 1;
      Lx = Ly = boost::lexical_cast<uint_t>(argv[1]);
      tmin = tmax = boost::lexical_cast<real_t>(argv[2]);
    } else if (argc == 4) {
      Lx = Ly = boost::lexical_cast<uint_t>(argv[1]);
      Jx = Jy = boost::lexical_cast<real_t>(argv[2]);
      tmin = tmax = boost::lexical_cast<real_t>(argv[3]);
    } else if (argc == 6) {
      Lx = boost::lexical_cast<uint_t>(argv[1]);
      Ly = boost::lexical_cast<uint_t>(argv[2]);
      Jx = boost::lexical_cast<real_t>(argv[3]);
      Jy = boost::lexical_cast<real_t>(argv[4]);
      tmin = tmax = boost::lexical_cast<real_t>(argv[5]);
    } else if (argc == 8) {
      Lx = boost::lexical_cast<uint_t>(argv[1]);
      Ly = boost::lexical_cast<uint_t>(argv[2]);
      Jx = boost::lexical_cast<real_t>(argv[3]);
      Jy = boost::lexical_cast<real_t>(argv[4]);
      tmin = boost::lexical_cast<real_t>(argv[5]);
      tmax = boost::lexical_cast<real_t>(argv[6]);
      dt = boost::lexical_cast<real_t>(argv[7]);
    } else throw(0);
  } catch (...) {
    std::cerr << "Usage: " << argv[0] << " L T\n";
    std::cerr << "       " << argv[0] << " L J T\n";
    std::cerr << "       " << argv[0] << " Lx Ly Jx Jy T\n";
    std::cerr << "       " << argv[0] << " Lx Ly Jx Jy Tmin Tmax dT\n";
    return 127;
  }
  std::cout << std::scientific << std::setprecision(std::numeric_limits<real_t>::digits10);
  for (real_t t = tmin; t < tmax + 0.01 * dt; t += dt) {
    auto result = ising::free_energy::square::finite_count(Lx, Ly, Jx, Jy, 1 / t);
    std::cout << Lx << ' ' << Ly << ' ' << Jx << ' ' << Jy << ' ' << t << ' '
              << std::get<0>(result) << ' ' << std::get<1>(result) << ' '
              << std::get<2>(result) << std::endl;
  }
}
