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

template <typename U>
inline typename U::root_type free_energy(U f, U) {
  return f.derivative(0);
}

template <typename U>
inline typename U::root_type energy(U f, U beta) {
  auto f0 = f.derivative(0);
  auto f1 = f.derivative(1);
  auto b0 = beta.derivative(0);
  return (f0 + b0 * f1);
}

template <typename U>
inline typename U::root_type specific_heat(U f, U beta) {
  auto f0 = f.derivative(0);
  auto f1 = f.derivative(1);
  auto f2 = f.derivative(2);
  auto b0 = beta.derivative(0);
  return -(b0 * b0 * (2 * f1 + b0 * f2));
}

template <typename T, typename U, typename W>
inline typename T::root_type free_energy(T f, U, W) {
  return f.derivative(0, 0);
}

template <typename T, typename U, typename W>
inline typename T::root_type energy(T f, U beta, W) {
  auto f00 = f.derivative(0, 0);
  auto f10 = f.derivative(1, 0);
  auto b0 = beta.derivative(0);
  return (f00 + b0 * f10);
}

template <typename T, typename U, typename W>
inline typename T::root_type specific_heat(T f, U beta, W) {
  auto f00 = f.derivative(0, 0);
  auto f10 = f.derivative(1, 0);
  auto f20 = f.derivative(2, 0);
  auto b0 = beta.derivative(0);
  return -(b0 * b0 * (2 * f10 + b0 * f20));
}

template <typename T, typename U, typename W>
inline typename T::root_type susceptibility(T f, U beta, W) {
  return -f.derivative(0, 2);
}

template <typename T, typename U, typename W>
inline typename T::root_type magnetization2(T f, U beta, W h) {
  auto b0 = beta.derivative(0);
  return susceptibility(f, beta, h) / b0;
}

}  // namespace free_energy
}  // namespace ising
