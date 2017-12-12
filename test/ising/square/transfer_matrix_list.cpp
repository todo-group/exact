/*****************************************************************************
*
* Copyright (C) 2011-2017 by Synge Todo <wistaria@phy.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Calculating free energy density of square lattice Ising model

#include <ising/square/transfer_matrix.hpp>
#include <lse/exp_number.hpp>
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
  std::vector<double> interaction_x(lat.get_length_x()), interaction_y(lat.get_length_x());
  int dim = 1 << lat.get_length_x();
  std::vector<double> v(dim);
  for (double t = t_min; t <= t_max; t += t_step) {
    double beta = 1 / t;
    lse::exp_double sum = 0;
    for (int i = 0; i < dim; ++i) {
      lse::exp_double weight = 1;
      for (int j = 0; j < dim; ++j) v[j] = 0;
      v[i] = 1;
      for (int y = 0; y < lat.get_length_y(); ++y) {
        for (int x = 0; x < lat.get_length_x(); ++x) {
          interaction_x[x] = interaction[2 * (Lx * y + x)];
          interaction_y[x] = interaction[2 * (Lx * y + x) + 1];
        }
        weight *= ising::square::transfer_matrix::product_D(beta, interaction_x, v);
        weight *= ising::square::transfer_matrix::product_U(beta, interaction_y, v);
      }
      weight *= v[i];
      sum += weight;
    }
    double f = -log(sum) / beta / lat.num_sites();
    std::cout << boost::format("%1% %2$.11e") % t % f << std::endl;
  }
}
