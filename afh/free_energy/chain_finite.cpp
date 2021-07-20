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
#include <Eigen/Dense>
#include "standards/exp_number.hpp"
#include "chain.hpp"

typedef std::size_t uint_t;
typedef Eigen::MatrixXd matrix_t;

struct options {
  uint_t L;
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

  // lattice
  std::vector<std::pair<uint_t, uint_t> > lattice;
  for (uint_t i = 0; i < opt.L; ++i) lattice.push_back(std::make_pair(i, (i+1) % opt.L));
  
  // generate Hamiltonian
  uint_t dim = 1 << opt.L;
  matrix_t hamiltonian(dim, dim);
  afh::free_energy::chain::generate(opt.L, lattice, hamiltonian);
  
  /* perform eigenvalue decomposition */
  Eigen::SelfAdjointEigenSolver<matrix_t> eigensolver(hamiltonian);
  auto evals = eigensolver.eigenvalues();
  std::cout << std::scientific << std::setprecision(11);
  std::cout << "# ground state energy/L: " << (evals(0) / opt.L) << std::endl;
  std::cout << "# first excitation gap: " << (evals(1) - evals(0)) << std::endl;

  std::cout << "# [T] [free energy/L] [energy/L] [entropy/L]\n";
  for (auto t = opt.Tmin; t < opt.Tmax + 1e-4 * opt.dT; t += opt.dT) {
    double beta = 1 / t;
    // calculate free energy and intternal energy
    standards::exp_double z = 0;
    standards::exp_double w = 0;
    for (uint_t i = 0; i < dim; ++i) {
      z += standards::exp_double::exp(-beta * evals(dim - i));
      w += evals(dim - i) * standards::exp_double::exp(-beta * evals(dim - i));
    }
    double f = - log(z) / (beta * opt.L);
    double e = w / z / opt.L;
    std::cout << t << ' ' << f << ' ' << e << ' ' << (beta * (e-f)) << std::endl;
  }
}
