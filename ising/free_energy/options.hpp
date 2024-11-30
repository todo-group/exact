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
  std::string Jx, Jy, Tmin, Tmax, dT;
  options2(unsigned argc, char *argv[])
      : valid(true), prec(15), Jx("1"), Jy("1") {
    if (argc == 1) {
      std::cerr << help(argv[0]);
      return;
    }
    for (unsigned i = 1; i < argc; ++i) {
      if (std::string(argv[i]) == "-h" || std::string(argv[i]) == "--help") {
        std::cerr << help(argv[0]);
        return;
      }
      switch (argv[i][0]) {
        case '-':
          switch (argv[i][1]) {
            case 'p':
              if (++i == argc) {
                std::cerr << help(argv[0]);
                return;
              }
              prec = atoi(argv[i]);
              break;
            default:
              std::cerr << help(argv[0]);
              return;
          }
          break;
        default:
          switch (argc - i) {
            case 1:
              Jx = Jy = "1";
              Tmin = Tmax = dT = argv[i];
              return;
            case 2:
              Jx = Jy = argv[i];
              Tmin = Tmax = dT = argv[i + 1];
              return;
            case 3:
              Jx = argv[i];
              Jy = argv[i + 1];
              Tmin = Tmax = dT = argv[i + 2];
              return;
            case 5:
              Jx = argv[i];
              Jy = argv[i + 1];
              Tmin = argv[i + 2];
              Tmax = argv[i + 3];
              dT = argv[i + 4];
              return;
            default:
              std::cerr << help(argv[0]);
              return;
          }
      }
    }
  }
  std::string help(char *prog) {
    valid = false;
    return std::string("Free energy of ferromagnetic Ising model\n") +
           "Usage: " + prog + " [-p prec] T\n" + "       " + prog +
           " [-p prec] J T\n" + "       " + prog + " [-p prec] Jx Jy T\n" +
           "       " + prog + " [-p prec] Jx Jy Tmin Tmax dT\n" +
           "Note: T can be specified as \"tc\" instead of real numbers\n";
  }
};

struct options2f {
  bool valid;
  unsigned int prec;
  unsigned long Lx, Ly;
  std::string Jx, Jy, Tmin, Tmax, dT;
  options2f(unsigned argc, char *argv[])
      : valid(true), prec(15), Jx("1"), Jy("1") {
    if (argc == 1) {
      std::cerr << help(argv[0]);
      return;
    }
    for (unsigned i = 1; i < argc; ++i) {
      if (std::string(argv[i]) == "-h" || std::string(argv[i]) == "--help") {
        std::cerr << help(argv[0]);
        return;
      }
      switch (argv[i][0]) {
        case '-':
          switch (argv[i][1]) {
            case 'p':
              if (++i == argc) {
                std::cerr << help(argv[0]);
                return;
              }
              prec = atoi(argv[i]);
              break;
            default:
              std::cerr << help(argv[0]);
              return;
          }
          break;
        default:
          switch (argc - i) {
            case 2:
              Lx = Ly = std::atol(argv[i]);
              Jx = Jy = "1";
              Tmin = Tmax = dT = argv[i + 1];
              return;
            case 3:
              Lx = Ly = std::atol(argv[i]);
              Jx = Jy = argv[i + 1];
              Tmin = Tmax = dT = argv[i + 2];
              return;
            case 5:
              Lx = std::atol(argv[i]);
              Ly = std::atol(argv[i + 1]);
              Jx = argv[i + 2];
              Jy = argv[i + 3];
              Tmin = Tmax = dT = argv[i + 4];
              return;
            case 7:
              Lx = std::atol(argv[i]);
              Ly = std::atol(argv[i + 1]);
              Jx = argv[i + 2];
              Jy = argv[i + 3];
              Tmin = argv[i + 4];
              Tmax = argv[i + 5];
              dT = argv[i + 6];
              return;
            default:
              std::cerr << help(argv[0]);
              return;
          }
      }
    }
  }
  std::string help(char *prog) {
    valid = false;
    return std::string("Free energy of ferromagnetic Ising model\n") +
           "Usage: " + prog + " [-p prec] L T\n" + "       " + prog +
           " [-p prec] L J T\n" + "       " + prog +
           " [-p prec] Lx Ly Jx Jy T\n" + "       " + prog +
           " [-p prec] Lx Ly Jx Jy Tmin Tmax dT\n" +
           "Note: T can be specified as \"tc\" instead of real numbers\n";
  }
};

struct options3 {
  bool valid;
  unsigned int prec;
  std::string Ja, Jb, Jc, Tmin, Tmax, dT;
  options3(unsigned argc, char *argv[])
      : valid(true), prec(15), Ja("1"), Jb("1"), Jc("1") {
    if (argc == 1) {
      std::cerr << help(argv[0]);
      return;
    }
    for (unsigned i = 1; i < argc; ++i) {
      if (std::string(argv[i]) == "-h" || std::string(argv[i]) == "--help") {
        std::cerr << help(argv[0]);
        return;
      }
      switch (argv[i][0]) {
        case '-':
          switch (argv[i][1]) {
            case 'p':
              if (++i == argc) {
                std::cerr << help(argv[0]);
                return;
              }
              prec = atoi(argv[i]);
              break;
            default:
              std::cerr << help(argv[0]);
              return;
          }
          break;
        default:
          switch (argc - i) {
            case 1:
              Ja = Jb = Jc = "1";
              Tmin = Tmax = dT = argv[i];
              return;
            case 2:
              Ja = Jb = Jc = argv[i];
              Tmin = Tmax = dT = argv[i + 1];
              return;
            case 4:
              Ja = argv[i];
              Jb = argv[i + 1];
              Jc = argv[i + 2];
              Tmin = Tmax = dT = argv[i + 3];
              return;
            case 6:
              Ja = argv[i];
              Jb = argv[i + 1];
              Jc = argv[i + 2];
              Tmin = argv[i + 3];
              Tmax = argv[i + 4];
              dT = argv[i + 5];
              return;
            default:
              std::cerr << help(argv[0]);
              return;
          }
      }
    }
  }
  std::string help(char *prog) {
    valid = false;
    return std::string("Free energy of ferromagnetic Ising model\n") +
           "Usage: " + prog + " [-p prec] T\n" + "       " + prog +
           " [-p prec] J T\n" + "       " + prog + " [-p prec] Ja Jb Jc T\n" +
           "       " + prog + " [-p prec] Ja Jb Jc Tmin Tmax dT\n" +
           "Note: T can be specified as \"tc\" instead of real numbers\n";
  }
};
