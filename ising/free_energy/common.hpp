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

#pragma once

namespace ising {
namespace free_energy {

template<typename U>
inline typename U::root_type free_energy(U f, U) {
  return f.derivative(0);
}

template<typename U>
inline typename U::root_type energy(U f, U beta) {
  return (f + beta * f.derivative(1)).derivative(0);
}

template<typename U>
inline typename U::root_type specific_heat(U f, U beta) {
  return -(beta * beta * (2 * f.derivative(1) + beta * f.derivative(2))).derivative(0);
}

}
}
