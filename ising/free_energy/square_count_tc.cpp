/*****************************************************************************
*
* Copyright (C) 2015-2021 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Free energy, energy, and specific heat of square lattice Ising model

#include <iomanip>
#include <iostream>
#include <boost/math/differentiation/autodiff.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "ising/mp_wrapper.hpp"
#include "square.hpp"
#include "ising/tc/square.hpp"

struct options {
  unsigned int prec;
  std::string Jx, Jy;
  unsigned long Lx, Ly;
  bool valid;
  options(unsigned int argc, char *argv[]) :
    prec(15), Jx("1"), Jy("1"), valid(true) {
    if (argc == 1) { valid = false; return; }
    for (unsigned i = 1; i < argc; ++i) {
      switch (argv[i][0]) {
      case '-' :
        switch (argv[i][1]) {
        case 'p' :
          if (++i == argc) { valid = false; return; }
          prec = std::atoi(argv[i]); break;
        default :
          valid = false; return;
        }
        break;
      default :
        switch (argc - i) {
        case 1:
          Lx = Ly = std::atol(argv[i]); return;
        case 2:
          Lx = Ly = std::atol(argv[i]);
          Jx = Jy = argv[i+1];
          return;
        case 4:
          Lx = std::atol(argv[i]);
          Ly = std::atol(argv[i+1]);
          Jx = argv[i+2];
          Jy = argv[i+3];
          return;
        default:
          valid = false; return;
        }
      }
    }
  }
};
  
template<typename T>
void calc(const options& opt) {
  typedef T real_t;
  namespace ifs = ising::free_energy::square;
  real_t Jx = convert<real_t>(opt.Jx);
  real_t Jy = convert<real_t>(opt.Jy);
  real_t tc = ising::tc::square(Jx, Jy);
  auto beta = boost::math::differentiation::make_fvar<real_t, 2>(1 / tc);
  auto f = ifs::finite_count(opt.Lx, opt.Ly, Jx, Jy, beta);
  std::cout << std::scientific << std::setprecision(std::numeric_limits<real_t>::digits10)
            << "# lattice: square\n"
            << "# precision: " << std::numeric_limits<real_t>::digits10 << std::endl
            << "# Lx Ly Jx Jy T 1/T F/N E/N C/N\n"
            << opt.Lx << ' ' << opt.Ly << ' ' << Jx << ' ' << Jy << ' '
            << tc << ' ' << (1 / tc) << ' '
            << ifs::free_energy(f, beta) << ' ' << ifs::energy(f, beta) << ' '
            << ifs::specific_heat(f, beta) << std::endl;
}

int main(int argc, char **argv) {
  namespace mp = boost::multiprecision;
  options opt(argc, argv);
  if (!opt.valid) {
    std::cerr << "Usage: " << argv[0] << " [-p prec] L\n"
              << "       " << argv[0] << " [-p prec] L J\n"
              << "       " << argv[0] << " [-p prec] Lx Ly Jx Jy\n"; return 127;
  }
  if (opt.prec <= std::numeric_limits<float>::digits10) {
    calc<float>(opt);
  } else if (opt.prec <= std::numeric_limits<double>::digits10) {
    calc<double>(opt);
  } else if (opt.prec <= std::numeric_limits<mp_wrapper<mp::cpp_dec_float_50>>::digits10) {
    calc<mp_wrapper<mp::cpp_dec_float_50>>(opt);
  } else if (opt.prec <= std::numeric_limits<mp_wrapper<mp::cpp_dec_float_100>>::digits10) {
    calc<mp_wrapper<mp::cpp_dec_float_100>>(opt);
  } else {
    std::cerr << "Error: Required precision is too high\n"; return 127;
  }
}
