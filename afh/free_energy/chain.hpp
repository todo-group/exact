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

// Generate Hamiltonian of antiferromagnetic Heisenberg Model

#pragma once

#include <vector>

namespace afh {
namespace free_energy {
namespace chain {
    
template<typename U, typename MATRIX>
void generate(U L, const std::vector<std::pair<U, U> >& lattice, MATRIX& mat) {
  typedef std::size_t uint_t;
  
  uint_t N = 1 << L;
  for (uint_t j = 0; j < N; ++j)
    for (uint_t i = 0; i < N; ++i)
      mat(i, j) = 0;

  for (uint_t l = 0; l < lattice.size(); ++l) {
    uint_t i = lattice[l].first;
    uint_t j = lattice[l].second;
    uint_t m1 = 1 << i;
    uint_t m2 = 1 << j;
    uint_t m3 = m1 + m2;
    for (uint_t k = 0; k < N; ++k) {
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

} // end namespace chain
} // end namespace free_energy
} // end namespace afh
