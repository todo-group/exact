/*****************************************************************************
*
* Copyright (C) 2011-2017 by Synge Todo <wistaria@phy.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Calculating free energy density of square lattice Ising model

#include <iomanip>
#include <iostream>
#include <string>
#include <lattice/graph.hpp>
#include "ising/counting.hpp"

int main(int argc, char **argv) {
  int L; // system size
  double J, H;
  double t_min, t_max, t_step;
  if (argc >=7) {
    L = std::stoi(argv[1]);
    J = std::stod(argv[2]);
    H = std::stod(argv[3]);
    t_min = std::stod(argv[4]);
    t_max = std::stod(argv[5]);
    t_step = std::stod(argv[6]);
  } else {
    std::cin >> L >> J >> H >> t_min >> t_max >> t_step;
  }
  std::cout << "# L = " << L << std::endl
            << "# J = " << J << std::endl
            << "# H = " << H << std::endl;
  lattice::graph lat = lattice::graph::simple(2, L);
  std::cout << std::scientific << std::setprecision(11);
  for (double t = t_min; t <= t_max; t += t_step) {
    double beta = 1 / t;
    double f = ising::counting::free_energy_density(beta, lat, J, H);
    std::cout << t << ' ' << f << std::endl;
  }
}
