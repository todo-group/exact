/*****************************************************************************
*
* Copyright (C) 2015-2021 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Critical temperature of triangular lattice Ising model

#include <iomanip>
#include <iostream>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "ising/mp_wrapper.hpp"
#include "triangular.hpp"

struct options {
  unsigned int prec;
  std::string Ja, Jb, Jc;
  bool valid;
  options(unsigned int argc, char *argv[]) :
    prec(15), Ja("1"), Jb("1"), Jc("1"), valid(true) {
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
          Ja = Jb = Jc = argv[i]; return;
        case 3:
          Ja = argv[i]; Jb = argv[i+1]; Jc = argv[i+2]; return;
        default:
          valid = false; return;
        }
      }
    }
  }
};

template<class T>
void calc(const std::string& Ja_in, const std::string& Jb_in, const std::string& Jc_in) {
  typedef T real_t;
  real_t Ja = convert<real_t>(Ja_in);
  real_t Jb = convert<real_t>(Jb_in);
  real_t Jc = convert<real_t>(Jc_in);
  auto tc = ising::tc::triangular(Ja, Jb, Jc);
  std::cout << std::scientific << std::setprecision(std::numeric_limits<real_t>::digits10)
            << "# lattice: triangular\n"
            << "# precision: " << std::numeric_limits<real_t>::digits10 << std::endl
            << "# Ja Jb Jc Tc 1/Tc" << std::endl
            << Ja << ' ' << Jb << ' ' << Jc << ' ' << tc << ' ' << (1 / tc) << std::endl;
}

int main(int argc, char **argv) {
  namespace mp = boost::multiprecision;
  options opt(argc, argv);
  if (!opt.valid) {
    std::cerr << "Usage: " << argv[0] << " [-p prec] [Ja [Jb Jc]]\n"; return 127;
  }
  if (opt.prec <= std::numeric_limits<float>::digits10) {
    calc<float>(opt.Ja, opt.Jb, opt.Jc);
  } else if (opt.prec <= std::numeric_limits<double>::digits10) {
    calc<double>(opt.Ja, opt.Jb, opt.Jc);
  } else if (opt.prec <= std::numeric_limits<mp_wrapper<mp::cpp_dec_float_50>>::digits10) {
    calc<mp_wrapper<mp::cpp_dec_float_50>>(opt.Ja, opt.Jb, opt.Jc);
  } else if (opt.prec <= std::numeric_limits<mp_wrapper<mp::cpp_dec_float_100>>::digits10) {
    calc<mp_wrapper<mp::cpp_dec_float_100>>(opt.Ja, opt.Jb, opt.Jc);
  } else {
    std::cerr << "Error: Required precision is too high\n"; return 127;
  }
}
