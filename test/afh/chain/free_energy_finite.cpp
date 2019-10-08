/*****************************************************************************
*
* Copyright (C) 2016-2018 by Synge Todo <wistaria@phy.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Calculating free energy density of quantum antiferomagnetic Heisenberg chain

#include <iostream>
#include <string>
#include <afh/chain/heisenberg.hpp>
#include <lse/exp_number.hpp>
#include <boost/format.hpp>
#include <boost/numeric/bindings/lapack/driver/syev.hpp>
#include <boost/numeric/bindings/lower.hpp>
#include <boost/numeric/bindings/ublas.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

using namespace boost::numeric;

typedef ublas::vector<double> vector_type;
typedef ublas::matrix<double, ublas::column_major> matrix_type;

int main(int argc, char** argv) {
  int L; // system size
  double t_min, t_max, t_step;
  if (argc >=5) {
    L = std::stoi(argv[1]);
    t_min = std::stod(argv[2]);
    t_max = std::stod(argv[3]);
    t_step = std::stod(argv[4]);
  } else {
    std::cin >> L >> t_min >> t_max >> t_step;
  }
  std::cout << "# L = " << L << std::endl;

  // lattice
  std::vector<std::pair<int, int> > lattice;
  for (int i = 0; i < L; ++i) lattice.push_back(std::make_pair(i, (i+1) % L));
  
  // generate Hamiltonian
  int dim = 1 << L;
  matrix_type hamiltonian(dim, dim);
  exact::afh::generate(L, lattice, hamiltonian);
  
  /* perform eigenvalue decomposition */
  vector_type evals(dim);
  char jobz = 'N'; /* Compute eigenvalues only */
  int info = bindings::lapack::syev(jobz, bindings::lower(hamiltonian), evals);
  if (info != 0) {
    std::cerr << "Error: LAPACK::dsyev failed\n";
    std::exit(1);
  }
  std::cout << boost::format("# ground state energy/L: %1$.11e") % (evals(0) / L) << std::endl;
  std::cout << boost::format("# first excitation gap: %1$.11e") % (evals(1) - evals(0)) << std::endl;

  std::cout << "# [T] [free energy/L] [energy/L] [entropy/L]\n";
  for (double t = t_min; t <= t_max; t += t_step) {
    double beta = 1 / t;
    // calculate free energy and internal energy
    lse::exp_double z = 0;
    lse::exp_double w = 0;
    for (int i = dim - 1; i >= 0; --i) {
      z += lse::exp_value(-beta * evals(i));
      w += evals(i) * lse::exp_value(-beta * evals(i));
    }
    double f = - log(z) / (beta * L);
    double e = w / z / L;
    std::cout << boost::format("%1% %2$.11e %3$.11e %4$.11e")
      % t % f % e % (beta * (e-f)) << std::endl;
  }
}
