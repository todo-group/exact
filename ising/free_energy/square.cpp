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
#include "options.hpp"
#include "square.hpp"

template<typename T>
void calc0(const options2& opt) {
  using namespace ising::free_energy;
  typedef T real_t;
  real_t Jx = convert<real_t>(opt.Jx);
  real_t Jy = convert<real_t>(opt.Jy);
  real_t Tmin = convert<real_t>(opt.Tmin);
  real_t Tmax = convert<real_t>(opt.Tmax);
  real_t dT = convert<real_t>(opt.dT);
  if (Tmin > Tmax) throw(std::invalid_argument("Tmax should be larger than Tmin"));
  if (dT <= 0) throw(std::invalid_argument("dT should be positive"));
  std::cout << std::scientific << std::setprecision(std::numeric_limits<real_t>::digits10)
            << "# lattice: square\n"
            << "# precision: " << std::numeric_limits<real_t>::digits10 << std::endl
            << "# Lx Ly Jx Jy T 1/T F/N E/N C/N\n";
  for (auto t = Tmin; t < Tmax + 1e-4 * dT; t += dT) {
    real_t beta = 1 / t;
    auto f = square::infinite(Jx, Jy, beta);
    std::cout << "inf inf " << Jx << ' ' << Jy << ' ' << t << ' ' << beta << ' '
              << f << " N/A N/A" << std::endl;
  }
}

template<typename T>
void calc(const options2& opt) {
  using namespace ising::free_energy;
  typedef T real_t;
  real_t Jx = convert<real_t>(opt.Jx);
  real_t Jy = convert<real_t>(opt.Jy);
  real_t Tmin = convert<real_t>(opt.Tmin);
  real_t Tmax = convert<real_t>(opt.Tmax);
  real_t dT = convert<real_t>(opt.dT);
  if (Tmin > Tmax) throw(std::invalid_argument("Tmax should be larger than Tmin"));
  if (dT <= 0) throw(std::invalid_argument("dT should be positive"));
  std::cout << std::scientific << std::setprecision(std::numeric_limits<real_t>::digits10)
            << "# lattice: square\n"
            << "# precision: " << std::numeric_limits<real_t>::digits10 << std::endl
            << "# Lx Ly Jx Jy T 1/T F/N E/N C/N\n";
  for (auto t = Tmin; t < Tmax + 1e-4 * dT; t += dT) {
    auto beta = boost::math::differentiation::make_fvar<real_t, 2>(1 / t);
    auto f = square::infinite(Jx, Jy, beta);
    std::cout << "inf inf " << Jx << ' ' << Jy << ' ' << t << ' ' << (1 / t) << ' '
              << free_energy(f, beta) << ' ' << energy(f, beta) << ' '
              << specific_heat(f, beta) << std::endl;
  }
}

int main(int argc, char **argv) {
  using namespace boost::multiprecision;
  options2 opt(argc, argv);
  if (!opt.valid) return 127;
  if (opt.prec <= std::numeric_limits<float>::digits10) {
    calc<float>(opt);
  } else if (opt.prec <= std::numeric_limits<double>::digits10) {
    calc<double>(opt);
  } else if (opt.prec <= std::numeric_limits<mp_wrapper<cpp_dec_float_50>>::digits10) {
    calc0<cpp_dec_float_50>(opt);
    // calc<mp_wrapper<cpp_dec_float_50>>(opt);
  } else if (opt.prec <= std::numeric_limits<mp_wrapper<cpp_dec_float_100>>::digits10) {
    calc0<cpp_dec_float_100>(opt);
    // calc<mp_wrapper<cpp_dec_float_100>>(opt);
  } else {
    std::cerr << "Error: Required precision is too high\n"; return 127;
  }
}
