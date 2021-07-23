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
#include <stdexcept>
#include <boost/math/differentiation/autodiff.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "ising/mp_wrapper.hpp"
#include "ising/tc/square.hpp"
#include "options.hpp"
#include "square.hpp"

template<typename T>
void calc(const options2f& opt) {
  using namespace ising::free_energy;
  typedef T real_t;
  real_t Jx = convert<real_t>(opt.Jx);
  real_t Jy = convert<real_t>(opt.Jy);
  real_t Tmin, Tmax, dT;
  if (opt.Tmin == "tc" || opt.Tmin == "Tc") {
    Tmin = Tmax = dT = ising::tc::square(Jx, Jy);
  } else {
    Tmin = convert<real_t>(opt.Tmin);
    Tmax = convert<real_t>(opt.Tmax);
    dT = convert<real_t>(opt.dT);
  }
  if (Tmin < 0 || Tmax < 0) throw(std::invalid_argument("Temperature should be positive"));
  if (Tmin > Tmax) throw(std::invalid_argument("Tmax should be larger than Tmin"));
  if (dT <= 0) throw(std::invalid_argument("dT should be positive"));
  std::cout << std::scientific << std::setprecision(std::numeric_limits<real_t>::digits10)
            << "# lattice: square\n"
            << "# precision: " << std::numeric_limits<real_t>::digits10 << std::endl
            << "# Lx Ly Jx Jy T 1/T F/N E/N C/N\n";
  for (auto t = Tmin; t < Tmax + 1e-4 * dT; t += dT) {
    auto beta = boost::math::differentiation::make_fvar<real_t, 2>(1 / t);
    auto f = square::finite_count(opt.Lx, opt.Ly, Jx, Jy, beta);
    std::cout << opt.Lx << ' ' << opt.Ly << ' ' << Jx << ' ' << Jy << ' '
              << t << ' ' << (1 / t) << ' '
              << free_energy(f, beta) << ' ' << energy(f, beta) << ' '
              << specific_heat(f, beta) << std::endl;
  }
}

int main(int argc, char **argv) {
  using namespace boost::multiprecision;
  options2f opt(argc, argv);
  if (!opt.valid) return 127;
  if (opt.prec <= std::numeric_limits<float>::digits10) {
    calc<float>(opt);
  } else if (opt.prec <= std::numeric_limits<double>::digits10) {
    calc<double>(opt);
  } else if (opt.prec <= std::numeric_limits<mp_wrapper<cpp_dec_float_50>>::digits10) {
    calc<mp_wrapper<cpp_dec_float_50>>(opt);
  } else if (opt.prec <= std::numeric_limits<mp_wrapper<cpp_dec_float_100>>::digits10) {
    calc<mp_wrapper<cpp_dec_float_100>>(opt);
  } else {
    std::cerr << "Error: Required precision is too high\n"; return 127;
  }
}
