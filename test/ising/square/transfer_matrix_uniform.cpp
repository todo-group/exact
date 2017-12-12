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
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>

int main(int argc, char **argv) {
  int L; // system size
  double J, H;
  double t_min, t_max, t_step;
  if (argc >=7) {
    L = boost::lexical_cast<int>(argv[1]);
    J = boost::lexical_cast<double>(argv[2]);
    H = boost::lexical_cast<double>(argv[3]);
    t_min = boost::lexical_cast<double>(argv[4]);
    t_max = boost::lexical_cast<double>(argv[5]);
    t_step = boost::lexical_cast<double>(argv[6]);
  } else {
    std::cin >> L >> J >> H >> t_min >> t_max >> t_step;
  }
  std::cout << "# L = " << L << std::endl
            << "# J = " << J << std::endl
            << "# H = " << H << std::endl;
  lattice::square lat(L, L);
  for (double t = t_min; t <= t_max; t += t_step) {
    double beta = 1 / t;
    double f = ising::square::transfer_matrix::free_energy_density(beta, lat, J, H);
    std::cout << boost::format("%1% %2$.11e") % t % f << std::endl;
  }
}
