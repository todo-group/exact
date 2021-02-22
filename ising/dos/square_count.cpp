/*****************************************************************************
*
* Copyright (C) 2015-2021 by Synge Todo <wistaria@phy.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Density of state of square lattice Ising model

#include<iostream>
#include<limits>
#include<vector>
#include <lattice/graph.hpp>

int main(int argc, char **argv) {
  typedef unsigned long uint_t;
  uint_t Lx, Ly;
  try {
    if (argc == 2) {
      Lx = Ly = std::stoi(argv[1]);
    } else if (argc == 3) {
      Lx = std::stoi(argv[1]);
      Ly = std::stoi(argv[2]);
    } else throw(0);
    if (Lx == 0 || Ly == 0) {
      std::cerr << "Error: Lx and Ly should be positive\n";
      throw(0);
    }
  } catch (...) {
    std::cerr << "Usage: " << argv[0] << " L\n";
    std::cerr << "       " << argv[0] << " Lx Ly\n";
    return 127;
  }

  auto basis = lattice::basis::simple(2);
  auto unitcell = lattice::unitcell(2);
  unitcell.add_site(lattice::coordinate(0, 0), 0);
  unitcell.add_bond(0, 0, lattice::offset(1, 0), 0);
  unitcell.add_bond(0, 0, lattice::offset(0, 1), 0);
  auto graph = lattice::graph(basis, unitcell, lattice::extent(Lx, Ly));
  if (graph.num_bonds() >= std::numeric_limits<uint_t>::digits) {
    std::cerr << "Error: system size is too large\n";
    return 127;
  }
  
  std::vector<uint_t> dos(graph.num_bonds() + 1, 0);
  uint_t num_states = 1 << graph.num_sites();
  for (uint_t c = 0; c < num_states; ++c) {
    uint_t energy = 0;
    for (uint_t b = 0; b < graph.num_bonds(); ++b) {
      uint_t ci = (c >> graph.source(b)) & 1;
      uint_t cj = (c >> graph.target(b)) & 1;
      energy += (ci ^ cj);
    }
    ++dos[energy];
  }
  for (uint_t i = 0; i <= graph.num_bonds(); i += 2) std::cout << dos[i] << ' ';
  std::cout << std::endl;
}
