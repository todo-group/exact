/*****************************************************************************
*
* Copyright (C) 2015-2020 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Critical temperature of square lattice Ising model

#include <iomanip>
#include <iostream>
#include <boost/lexical_cast.hpp>
#include "square.hpp"

int main(int argc, char **argv) {
  typedef double real_t;
  real_t Jx, Jy;
  try {
    if (argc == 1) {
      Jx = Jy = 1;
    } else if (argc == 2) {
      Jx = Jy = boost::lexical_cast<real_t>(argv[1]);
    } else if (argc == 3) {
      Jx = boost::lexical_cast<real_t>(argv[1]);
      Jy = boost::lexical_cast<real_t>(argv[2]);
    } else throw(0);
  } catch (...) {
    std::cerr << "Usage: " << argv[0] << " [Jx [Jy]]\n";
    return 127;
  }
  std::cout << std::scientific << std::setprecision(std::numeric_limits<real_t>::digits10)
            << ising::tc::square(Jx, Jy) << std::endl;
}
