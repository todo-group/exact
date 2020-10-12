/*****************************************************************************
*
* Copyright (C) 2015-2020 by Synge Todo <wistaria@phy.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Calculating free energy, energy, and specific heat of square lattice Ising model

#include <iomanip>
#include <iostream>
#include <boost/lexical_cast.hpp>
#include "square/infinite.hpp"

int main(int argc, char **argv) {
  typedef double real_t;
  real_t Jx, Jy, t_min, t_max, t_step;
  if (argc == 6) {
    Jx = boost::lexical_cast<real_t>(argv[1]);
    Jy = boost::lexical_cast<real_t>(argv[2]);
    t_min = boost::lexical_cast<real_t>(argv[3]);
    t_max = boost::lexical_cast<real_t>(argv[4]);
    t_step = boost::lexical_cast<real_t>(argv[5]);
  } else if (argc == 1) {
    std::cin >> Jx >> Jy >> t_min >> t_max >> t_step;
  } else {
    std::cerr << "Usage: " << argv[0] << " [Jx Jy t_min t_max t_step]\n";
    return 127;
  }
  std::cout << std::scientific << std::setprecision(std::numeric_limits<real_t>::digits10);
  std::cout << "# square lattice Ising model\n";
  std::cout << "# Jx, Jy, T, free energy density, energy density, specific heat\n";
  for (real_t t = t_min; t <= t_max; t += t_step) {
    real_t beta = 1 / t;
    auto result = ising::square::infinite(beta, Jx, Jy);
    std::cout << Jx << ' ' << Jy << ' ' << t << ' ' << std::get<0>(result) << ' '
              << std::get<1>(result) << ' ' << std::get<2>(result) << std::endl;
  }
}
