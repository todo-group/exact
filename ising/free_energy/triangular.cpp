/*****************************************************************************
*
* Copyright (C) 2015-2020 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Free energy, energy, and specific heat of triangular lattice Ising model

#include <iomanip>
#include <iostream>
#include <boost/lexical_cast.hpp>
#include "triangular.hpp"

int main(int argc, char **argv) {
  typedef double real_t;
  real_t Ja, Jb, Jc, t;
  try {
    if (argc == 2) {
      Ja = Jb = Jc = 1;
      t = boost::lexical_cast<real_t>(argv[1]);
    } else if (argc == 3) {
      Ja = Jb = Jc = boost::lexical_cast<real_t>(argv[1]);
      t = boost::lexical_cast<real_t>(argv[2]);
    } else if (argc == 5) {
      Ja = boost::lexical_cast<real_t>(argv[1]);
      Jb = boost::lexical_cast<real_t>(argv[2]);
      Jc = boost::lexical_cast<real_t>(argv[3]);
      t = boost::lexical_cast<real_t>(argv[4]);
    } else throw(0);
  } catch (...) {
    std::cerr << "Usage: " << argv[0] << " T\n";
    std::cerr << "       " << argv[0] << " J T\n";
    std::cerr << "       " << argv[0] << " Ja Jb Jc T\n";
    return 127;
  }
  std::cout << std::scientific << std::setprecision(std::numeric_limits<real_t>::digits10);
  auto result = ising::free_energy::triangular::infinite(1 / t, Ja, Jb, Jc);
  std::cout << Ja << ' ' << Jb << ' ' << Jc << ' ' << t << ' ' << std::get<0>(result) << ' '
            << std::get<1>(result) << ' ' << std::get<2>(result) << std::endl;
}
