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

#pragma once
#include <iostream>
#include <string>

struct options {
  bool valid;
  unsigned Lx, Ly;
  options(unsigned argc, char *argv[]) : valid(true) {
    switch (argc) {
    case 2:
      if (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help")
        std::cerr << help(argv[0]);
      else
        Lx = Ly = std::atol(argv[1]);
      return;
    case 3:
      Lx = std::atol(argv[1]);
      Ly = std::atol(argv[2]);
      return;
    default:
      std::cerr << help(argv[0]);
      return;
    }
  }
  std::string help(char *prog) {
    valid = false;
    return std::string("Density of state of square lattice Ising model\n") +
      "Usage: " + prog + " L\n" +
      "       " + prog + " Lx Ly\n";
  }
};
