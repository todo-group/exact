/*****************************************************************************
*
* Copyright (C) 2015-2020 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Critical temperature of hexagonal lattice Ising model

#include <iomanip>
#include <iostream>
#include <boost/lexical_cast.hpp>
#include "hexagonal.hpp"

int main(int argc, char **argv) {
  typedef double real_t;
  real_t Ja, Jb, Jc;
  try {
    if (argc == 1) {
      Ja = Jb = Jc = 1;
    } else if (argc == 2) {
      Ja = Jb = Jc = boost::lexical_cast<real_t>(argv[1]);
    } else if (argc == 4) {
      Ja = boost::lexical_cast<real_t>(argv[1]);
      Jb = boost::lexical_cast<real_t>(argv[2]);
      Jc = boost::lexical_cast<real_t>(argv[3]);
    } else throw(0);
  } catch (...) {
    std::cerr << "Usage: " << argv[0] << " [Ja [Jb Jc]]\n";
    return 127;
  }
  std::cout << std::scientific << std::setprecision(std::numeric_limits<real_t>::digits10)
            << Ja << ' ' << Jb << ' ' << Jc << ' '
            << ising::tc::hexagonal(Ja, Jb, Jc) << std::endl;
}
