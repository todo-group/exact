/*****************************************************************************
*
* Copyright (C) 2016 by Synge Todo <wistaria@phy.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Generate Hamiltonian of antiferromagnetic Heisenberg Model

#include <vector>

template<typename MATRIX>
void generate(int L, const std::vector<std::pair<int, int> >& lattice,
              MATRIX& mat) {
  int N = 1 << L;
  mat.clear();

  for (int l = 0; l < lattice.size(); ++l) {
    int i = lattice[l].first;
    int j = lattice[l].second;
    int m1 = 1 << i;
    int m2 = 1 << j;
    int m3 = m1 + m2;
    for (int k=0; k<N; ++k) {
      if (((k & m3) == m1) || ((k & m3) == m2)) {
        // when (bit i == 1, bit j == 0) or (bit i == 0, bit j == 1)
        mat(k^m3, k) += 0.5;
        mat(k, k) -= 0.25;
      } else {
        mat(k, k) += 0.25;
      }
    }
  }
}
