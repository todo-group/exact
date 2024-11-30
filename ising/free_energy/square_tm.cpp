/*
   Copyright (C) 2015-2021 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

// Free energy, energy, and specific heat of square lattice Ising model

#include <iomanip>
#include <iostream>
#include <stdexcept>
#include "ising/tc/square.hpp"
#include "ising/mp_wrapper.hpp"
#include "common.hpp"
#include "options.hpp"
#include "square_tm.hpp"

using namespace ising::free_energy;

template<typename T>
void calc(const options2f& opt) {
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
            << "# Lx Ly Jx Jy T 1/T F/N E/N C/N M2/N2\n";
  for (auto t = Tmin; t < Tmax + 1e-4 * dT; t += dT) {
    real_t h = 0;
    auto f = square::transfer_matrix<real_t>::calc(opt.Lx, opt.Ly, Jx, Jy, 1/t, h);
    typename square::transfer_matrix<real_t>::result_t beta;
    beta.set(0, 1 / t);
    std::cout << opt.Lx << ' ' << opt.Ly << ' ' << Jx << ' ' << Jy << ' '
              << t << ' ' << (1 / t) << ' '
              << free_energy(f, beta, h) << ' ' << energy(f, beta, h) << ' '
              << specific_heat(f, beta, h) << ' ' << magnetization2(f, beta, h)
              << std::endl;
  }
}

int main(int argc, char **argv) {
  options2f opt(argc, argv);
  if (!opt.valid) return 127;
  if (opt.prec <= std::numeric_limits<float>::digits10) {
    calc<float>(opt);
  } else if (opt.prec <= std::numeric_limits<double>::digits10) {
    calc<double>(opt);
  } else {
    std::cerr << "Error: Required precision is too high\n"; return 127;
  }
}
