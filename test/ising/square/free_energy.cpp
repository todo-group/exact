/*****************************************************************************
*
* Copyright (C) 2015-2016 by Synge Todo <wistaria@phy.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Calculating free energy density of square lattice Ising model

#include <exact/ising/square/free_energy.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>

int main(int argc, char **argv) {
  double t_min, t_max, t_step;
  int Nint;
  if (argc >=5) {
    t_min = boost::lexical_cast<double>(argv[1]);
    t_max = boost::lexical_cast<double>(argv[2]);
    t_step = boost::lexical_cast<double>(argv[3]);
    Nint = boost::lexical_cast<double>(argv[4]);
  } else {
    std::cin >> t_min >> t_max >> t_step >> Nint;
  }
  std::cout << "# Nint = " << Nint << std::endl;
  for (double t = t_min; t <= t_max; t += t_step) {
    double beta = 1 / t;
    double f = ising::square::free_energy_density(beta, 1, 1, Nint);
    std::cout << boost::format("%1% %2$.11e") % t % f << std::endl;
  }
}
