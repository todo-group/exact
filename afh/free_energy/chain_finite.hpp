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

#pragma once

#include <Eigen/Dense>
#include "standards/exp_number.hpp"
#include "chain.hpp"

namespace afh { namespace free_energy { namespace chain {
    
class finite {
  typedef std::size_t uint_t;
  typedef Eigen::VectorXd vector_t;
  typedef Eigen::MatrixXd matrix_t;
public:
  finite(uint_t L) : L_(L), eigenvalues_(1 << L_) {
    // lattice
    std::vector<std::pair<uint_t, uint_t> > lattice;
    for (uint_t i = 0; i < L_; ++i) lattice.push_back(std::make_pair(i, (i+1) % L_));
    
    // generate Hamiltonian
    uint_t dim = 1 << L_;
    matrix_t hamiltonian(dim, dim);
    afh::free_energy::chain::generate(L_, lattice, hamiltonian);
  
    /* perform eigenvalue decomposition */
    Eigen::SelfAdjointEigenSolver<matrix_t> eigensolver(hamiltonian);
    eigenvalues_ = eigensolver.eigenvalues();
  }
  double gs_energy() const { return eigenvalues_(0) / L_; }
  double gap() const { return eigenvalues_(1) - eigenvalues_(0); }
  std::tuple<double, double, double> free_energy(double t) const {
    // calculate free energy and intternal energy
    uint_t dim = 1 << L_;
    double beta = 1 / t;
    standards::exp_double z = 0;
    standards::exp_double w = 0;
    for (uint_t i = 0; i < dim; ++i) {
      uint_t j = dim - i - 1;
      z += standards::exp_double::exp(-beta * eigenvalues_(j));
      w += eigenvalues_(j) * standards::exp_double::exp(-beta * eigenvalues_(j));
    }
    double f = - log(z) / (beta * L_);
    double e = w / z / L_;
    return std::make_tuple(f, e, beta * (e-f));
  }
private:
  uint_t L_;
  vector_t eigenvalues_;
};

} } }
