/*
   Copyright (C) 2016-2021 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
                              Chihiro Kondo <chihiro.kondo@phys.s.u-tokyo.ac.jp>

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

// Density of state of square lattice Ising model

#include <iostream>
#include "options.hpp"
#include "square.hpp"

int main(int argc, char **argv) {
  options opt(argc, argv);
  if (!opt.valid) return 127;
  if (2 * opt.Lx * opt.Ly <= 128) {
    auto dos = ising::dos::square::finite<128, 100>(opt.Lx, opt.Ly);
    for (auto v : dos) std::cout << v << ' '; std::cout << std::endl;
  } else if (2 * opt.Lx * opt.Ly <= 512) {
    auto dos = ising::dos::square::finite<512, 100>(opt.Lx, opt.Ly);
  for (auto v : dos) std::cout << v << ' '; std::cout << std::endl;
  } else if (2 * opt.Lx * opt.Ly <= 2048) {
    auto dos = ising::dos::square::finite<2048, 100>(opt.Lx, opt.Ly);
    for (auto v : dos) std::cout << v << ' '; std::cout << std::endl;
  } else if (2 * opt.Lx * opt.Ly <= 4608) {
    auto dos = ising::dos::square::finite<4608, 100>(opt.Lx, opt.Ly);
    for (auto v : dos) std::cout << v << ' '; std::cout << std::endl;
  } else {
    std::cerr << "Error: Lx * Ly is too large\n";
    return 127;
  }
}
