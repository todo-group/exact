/*****************************************************************************
*
* Copyright (C) 2015-2021 by Synge Todo <wistaria@phy.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

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
