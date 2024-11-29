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
#include "ising/square/transfer_matrix.hpp"

int main() {
  int Lx, Ly; // system size
  std::vector<double> inter, field;
  double t_min, t_max, t_step;

  std::cin >> Lx >> Ly;
  inter.resize(2 * Lx * Ly);
  for (std::size_t b = 0; b < 2 * Lx * Ly; ++b) std::cin >> inter[b];
  field.resize(Lx * Ly);
  for (std::size_t s = 0; s < Lx * Ly; ++s) std::cin >> field[s];
  std::cin >> t_min >> t_max >> t_step;

  std::cout << "# Lx = " << Lx <<std::endl << "# Ly = " << Ly <<std::endl;
  std::cout << std::scientific << std::setprecision(11);
  for (double t = t_min; t <= t_max; t += t_step) {
    double beta = 1 / t;
    double f = ising::square::transfer_matrix::free_energy_density(beta, Lx, Ly, inter, field);
    std::cout << t << ' ' << f << std::endl;
  }
}
