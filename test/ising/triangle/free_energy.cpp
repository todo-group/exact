/*****************************************************************************
*
* Copyright (C) 2015-2020 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Calculating free energy, energy, and specific heat of triangular lattice Ising model

#include <iomanip>
#include <iostream>
#include <boost/lexical_cast.hpp>
#include "ising/free_energy/triangular.hpp"

int main(int argc, char **argv) {
  typedef double real_t;
  real_t Ja, Jb, Jc, t_min, t_max, t_step;
  if (argc == 7) {
    Ja = boost::lexical_cast<real_t>(argv[1]);
    Jb = boost::lexical_cast<real_t>(argv[2]);
    Jc = boost::lexical_cast<real_t>(argv[3]);
    t_min = boost::lexical_cast<real_t>(argv[4]);
    t_max = boost::lexical_cast<real_t>(argv[5]);
    t_step = boost::lexical_cast<real_t>(argv[6]);
  } else if (argc == 1) {
    std::cin >> Ja >> Jb >> Jc >> t_min >> t_max >> t_step;
  } else {
    std::cerr << "Usage: " << argv[0] << " [Ja Jb Jc t_min t_max t_step]\n";
    return 127;
  }
  std::cout << std::scientific << std::setprecision(std::numeric_limits<real_t>::digits10);
  std::cout << "# triangular lattice Ising model\n";
  std::cout << "# Ja, Jb, Jc, T, free energy density, energy density, specific heat\n";
  for (real_t t = t_min; t <= t_max; t += t_step) {
    auto beta = boost::math::differentiation::make_fvar<real_t, 2>(1 / t);
    auto result = ising::free_energy::triangular::infinite(Ja, Jb, Jc, beta);
    std::cout << Ja << ' ' << Jb << ' ' << Jc << ' ' << t << ' ' << result << std::endl;
  }
}
