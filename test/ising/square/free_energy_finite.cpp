/*****************************************************************************
*
* Copyright (C) 2011-2016 by Synge Todo <wistaria@phy.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Calculating free energy density of square lattice Ising model

#include <iomanip>
#include <iostream>
#include <string>
#include "ising/square/free_energy_finite.hpp"

int main(int argc, char **argv) {
  int L; // system size
  double t_min, t_max, t_step;
  if (argc >=5) {
    L = std::stoi(argv[1]);
    t_min = std::stod(argv[2]);
    t_max = std::stod(argv[3]);
    t_step = std::stod(argv[4]);
  } else {
    std::cin >> L >> t_min >> t_max >> t_step;
  }
  std::cout << "# L = " << L << std::endl;
  std::cout << std::scientific << std::setprecision(11);
  for (double t = t_min; t <= t_max; t += t_step) {
    double beta = 1 / t;
    double f = ising::square::free_energy_density(beta, 1, 1, L, L);
    std::cout << t << ' ' << f << std::endl;
  }
}
