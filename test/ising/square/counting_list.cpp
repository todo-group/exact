/*****************************************************************************
*
* Copyright (C) 2011-2017 by Synge Todo <wistaria@phy.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Calculating free energy density of square lattice Ising model

#include <ising/square/counting.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>

int main(int argc, char **argv) {
  int Lx, Ly; // system size
  std::vector<double> interaction;
  double t_min, t_max, t_step;

  std::cin >> Lx >> Ly;
  lattice::square lat(Lx, Ly);
  interaction.resize(lat.num_bonds());
  for (int b = 0; b < lat.num_bonds(); ++b) {
    std::cin >> interaction[b];
  }
  std::cin >> t_min >> t_max >> t_step;
  
  std::cout << "# Lx = " << Lx <<std::endl << "# Ly = " << Ly <<std::endl;
  for (double t = t_min; t <= t_max; t += t_step) {
    double beta = 1 / t;
    double f = ising::square::counting::free_energy_density(beta, lat, interaction);
    std::cout << boost::format("%1% %2$.11e") % t % f << std::endl;
  }
}