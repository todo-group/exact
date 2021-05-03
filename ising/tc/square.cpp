/*****************************************************************************
*
* Copyright (C) 2015-2021 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Critical temperature of square lattice Ising model

#include <iomanip>
#include <iostream>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "ising/mp_wrapper.hpp"
#include "square.hpp"

struct options {
  unsigned int prec;
  std::string Jx, Jy;
  bool valid;
  options(unsigned int argc, char *argv[]) :
    prec(15), Jx("1"), Jy("1"), valid(true) {
    for (unsigned i = 1; i < argc; ++i) {
      switch (argv[i][0]) {
      case '-' :
        switch (argv[i][1]) {
        case 'p' :
          if (++i == argc) { valid = false; return; }
          prec = atoi(argv[i]); break;
        default :
          valid = false; return;
        }
        break;
      default :
        switch (argc - i) {
        case 1:
          Jx = Jy = argv[i]; return;
        case 2:
          Jx = argv[i]; Jy = argv[i+1]; return;
        default:
          valid = false; return;
        }
      }
    }
  }
};
  
template<class T>
void calc(const std::string& Jx_in, const std::string& Jy_in) {
  typedef T real_t;
  real_t Jx = convert<real_t>(Jx_in);
  real_t Jy = convert<real_t>(Jy_in);
  auto tc = ising::tc::square(Jx, Jy);
  std::cout << std::scientific << std::setprecision(std::numeric_limits<real_t>::digits10)
            << "# lattice: square\n"
            << "# precision: " << std::numeric_limits<real_t>::digits10 << std::endl
            << "# Jx Jy Tc 1/Tc" << std::endl
            << Jx << ' ' << Jy << ' ' << tc << ' ' << (1 / tc) << std::endl;
}

int main(int argc, char **argv) {
  namespace mp = boost::multiprecision;
  options opt(argc, argv);
  if (!opt.valid) {
    std::cerr << "Usage: " << argv[0] << " [-p prec] [Jx [Jy]]\n"; return 127;
  }
  if (opt.prec <= std::numeric_limits<float>::digits10) {
    calc<float>(opt.Jx, opt.Jy);
  } else if (opt.prec <= std::numeric_limits<double>::digits10) {
    calc<double>(opt.Jx, opt.Jy);
  } else if (opt.prec <= std::numeric_limits<mp_wrapper<mp::cpp_dec_float_50>>::digits10) {
    calc<mp_wrapper<mp::cpp_dec_float_50>>(opt.Jx, opt.Jy);
  } else if (opt.prec <= std::numeric_limits<mp_wrapper<mp::cpp_dec_float_100>>::digits10) {
    calc<mp_wrapper<mp::cpp_dec_float_100>>(opt.Jx, opt.Jy);
  } else {
    std::cerr << "Error: Required precision is too high\n"; return 127;
  }
}
