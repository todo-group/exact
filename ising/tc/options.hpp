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

struct options2 {
  bool valid;
  unsigned int prec;
  std::string Jx, Jy;
  options2(unsigned argc, char *argv[]) : valid(true), prec(15), Jx("1"), Jy("1") {
    for (unsigned i = 1; i < argc; ++i) {
      if (argv[i] == "-h" || argv[i] == "--help") { std::cerr << help(argv[0]); return; }
      switch (argv[i][0]) {
      case '-' :
        switch (argv[i][1]) {
        case 'p' :
          if (++i == argc) { std::cerr << help(argv[0]); return; }
          prec = atoi(argv[i]); break;
        default :
          std::cerr << help(argv[0]); return;
        }
        break;
      default :
        switch (argc - i) {
        case 1:
          Jx = Jy = argv[i]; return;
        case 2:
          Jx = argv[i]; Jy = argv[i+1]; return;
        default:
          std::cerr << help(argv[0]); return;
        }
      }
    }
  }
  std::string help(char *prog) {
    valid = false;
    return std::string("Critical temperature of ferromagnetic Ising model\n") +
      "Usage: " + prog + " [-p prec] [Jx [Jy]]\n";
  }
};

struct options3 {
  bool valid;
  unsigned int prec;
  std::string Ja, Jb, Jc;
  options3(unsigned argc, char *argv[]) : valid(true), prec(15), Ja("1"), Jb("1"), Jc("1") {
    for (unsigned i = 1; i < argc; ++i) {
      if (argv[i] == "-h" || argv[i] == "--help") { std::cerr << help(argv[0]); return; }
      switch (argv[i][0]) {
      case '-' :
        switch (argv[i][1]) {
        case 'p' :
          if (++i == argc) { std::cerr << help(argv[0]); return; }
          prec = atoi(argv[i]); break;
        default :
          std::cerr << help(argv[0]); return;
        }
        break;
      default :
        switch (argc - i) {
        case 1:
          Ja = Jb = argv[i]; return;
        case 3:
          Ja = argv[i]; Jb = argv[i+1]; Jc = argv[i+2]; return;
        default:
          std::cerr << help(argv[0]); return;
        }
      }
    }
  }
  std::string help(char *prog) {
    valid = false;
    return std::string("Critical temperature of ferromagnetic Ising model\n") +
      "Usage: " + prog + " [-p prec] [Ja [Jb Jc]]\n";
  }
};
