/*****************************************************************************
*
* Copyright (C) 2011-2017 by Synge Todo <wistaria@phy.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Calculating free energy density of square lattice Ising model

#include <iomanip>
#include <iostream>
#include <string>
#include "ising/square/transfer_matrix.hpp"

int main(int argc, char **argv) {
  int L; // system size
  double t;
  if (argc >=3) {
    L = std::stoi(argv[1]);
    t = std::stod(argv[2]);
  } else {
    std::cin >> L >> t;
  }
  std::cout << "# L = " << L << std::endl;
  std::vector<double> inter(L, 1.0);
  std::vector<double> field(L, 0.1);
  int dim = 1 << L;
  std::vector<double> v(dim);
  double beta = 1 / t;
  std::cout << std::scientific << std::setprecision(5);
  
  std::cout << "# matrix D:\n";
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) v[j] = 0;
    v[i] = 1;
    ising::square::transfer_matrix::exp_double weight =
      ising::square::transfer_matrix::product_D(beta, inter, field, v);
    for (int j = 0; j < dim; ++j)
      std::cout << ' ' << double(weight * v[j]);
    std::cout << std::endl;
  }

  std::cout << "# matrix U:\n";
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) v[j] = 0;
    v[i] = 1;
    standards::exp_number<double> weight = ising::square::transfer_matrix::product_U(beta, inter, v);
    for (int j = 0; j < dim; ++j)
      std::cout << ' ' << double(weight * v[j]);
    std::cout << std::endl;
  }

  std::cout << "# matrix T = D^{1/2} U D^{1/2}:\n";
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) v[j] = 0;
    v[i] = 1;
    standards::exp_number<double> weight = 1;
    weight *= ising::square::transfer_matrix::product_D(beta / 2, inter, field, v);
    weight *= ising::square::transfer_matrix::product_U(beta, inter, v);
    weight *= ising::square::transfer_matrix::product_D(beta / 2, inter, field, v);
    for (int j = 0; j < dim; ++j)
      std::cout << ' ' << double(weight * v[j]);
    std::cout << std::endl;
  }
}
