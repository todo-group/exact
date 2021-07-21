/*
   Copyright (C) 2016-2021 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>

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

// Calculating free energy density of quantum antiferomagnetic Heisenberg chain

#include <iomanip>
#include <iostream>
#include <string>
#include "chain_finite.hpp"

struct options {
  unsigned L;
  double Tmin, Tmax, dT;
  bool valid;
  options(unsigned int argc, char *argv[]) : valid(true) {
    if (argc == 1) { valid = false; return; }
    switch (argc) {
    case 3:
      L = std::atoi(argv[1]);
      Tmin = Tmax = dT = std::atof(argv[2]);
      return;
    case 5:
      L = std::atoi(argv[1]);
      Tmin = std::atof(argv[2]);
      Tmax = std::atof(argv[3]);
      dT = std::atof(argv[4]);
      return;
    default:
      valid = false; return;
    }
  }
};

int main(int argc, char** argv) {
  options opt(argc, argv);
  if (!opt.valid) {
    std::cerr << "Usage: " << argv[0] << " L T\n"
              << "       " << argv[0] << " L Tmin Tmax dT\n";
    return 127;
  }
  if (opt.Tmin > opt.Tmax) throw(std::invalid_argument("Tmax should be larger than Tmin"));
  if (opt.dT <= 0) throw(std::invalid_argument("dT should be positive"));
  std::cout << "# L = " << opt.L << std::endl;

  afh::free_energy::chain::finite solver(opt.L);
  std::cout << std::scientific << std::setprecision(11);
  std::cout << "# ground state energy/L: " << solver.gs_energy() << std::endl;
  std::cout << "# first excitation gap: " << solver.gap() << std::endl;

  std::cout << "# [T] [free energy/L] [energy/L] [entropy/L]\n";
  for (auto t = opt.Tmin; t < opt.Tmax + 1e-4 * opt.dT; t += opt.dT) {
    double f, e, s;
    std::tie(f, e, s) = solver.free_energy(t);
    std::cout << t << ' ' << f << ' ' << e << ' ' << s << std::endl;
  }
}
