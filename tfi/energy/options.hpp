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
  unsigned int prec;
  std::string J, GammaMin, GammaMax, dGamma;
  options(unsigned argc, char *argv[]) : valid(true), prec(15), J("1"), dGamma("1") {
    if (argc == 1) { std::cerr << help(argv[0]); return; }
    for (unsigned i = 1; i < argc; ++i) {
      if (std::string(argv[i]) == "-h" || std::string(argv[i]) == "--help") { std::cerr << help(argv[0]); return; }
      if (std::string(argv[i]) == "-p") {
        if (++i == argc) { std::cerr << help(argv[0]); return; }
        prec = atoi(argv[i]);
      } else {
        switch (argc - i) {
        case 1:
          GammaMin = GammaMax = argv[i];
          return;
        case 2:
          J = argv[i];
          GammaMin = GammaMax = argv[i+1];
          return;
        case 4:
          J = argv[i];
          GammaMin = argv[i+1];
          GammaMax = argv[i+2];
          dGamma = argv[i+3];
          return;
        default:
          std::cerr << help(argv[0]); return;
        }
      }
    }
  }
  std::string help(char *prog) {
    valid = false;
    return std::string("Free energy of ferromagnetic Ising model\n") +
      "Usage: " + prog + " [-p prec] Gamma\n" +
      "       " + prog + " [-p prec] J Gamma\n" +
      "       " + prog + " [-p prec] J GammaMin GammaMax dGamma\n" +
      "Note: Gamma can be specified as \"GammaC\" instead of real numbers\n";
  }
};
