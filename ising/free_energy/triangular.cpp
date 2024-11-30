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

// Free energy, energy, and specific heat of triangular lattice Ising model

#include <iomanip>
#include <iostream>
#include <boost/math/differentiation/autodiff.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "ising/mp_wrapper.hpp"
#include "triangular.hpp"

struct options {
  unsigned int prec;
  std::string Ja, Jb, Jc, Tmin, Tmax, dT;
  bool valid;
  options(unsigned int argc, char *argv[])
      : prec(15), Ja("1"), Jb("1"), Jc("1"), valid(true) {
    if (argc == 1) {
      valid = false;
      return;
    }
    for (unsigned i = 1; i < argc; ++i) {
      switch (argv[i][0]) {
        case '-':
          switch (argv[i][1]) {
            case 'p':
              if (++i == argc) {
                valid = false;
                return;
              }
              prec = std::atoi(argv[i]);
              break;
            default:
              valid = false;
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
              valid = false;
              return;
          }
      }
    }
  }
};

template <typename T>
void calc(const options &opt) {
  using namespace ising::free_energy;
  typedef T real_t;
  real_t Ja = convert<real_t>(opt.Ja);
  real_t Jb = convert<real_t>(opt.Jb);
  real_t Jc = convert<real_t>(opt.Jc);
  real_t Tmin = convert<real_t>(opt.Tmin);
  real_t Tmax = convert<real_t>(opt.Tmax);
  real_t dT = convert<real_t>(opt.dT);
  if (Tmin > Tmax)
    throw(std::invalid_argument("Tmax should be larger than Tmin"));
  if (dT <= 0)
    throw(std::invalid_argument("dT should be positive"));
  std::cout << std::scientific
            << std::setprecision(std::numeric_limits<real_t>::digits10)
            << "# lattice: triangular\n"
            << "# precision: " << std::numeric_limits<real_t>::digits10
            << std::endl
            << "# Lx Ly Ja Jb Jc T 1/T F/N E/N C/N\n";
  for (auto t = Tmin; t < Tmax + 1e-4 * dT; t += dT) {
    auto beta = boost::math::differentiation::make_fvar<real_t, 2>(1 / t);
    auto f = triangular::infinite(Ja, Jb, Jc, beta);
    std::cout << "inf inf " << Ja << ' ' << Jb << ' ' << Jc << ' ' << t << ' '
              << (1 / t) << ' ' << free_energy(f, beta) << ' '
              << energy(f, beta) << ' ' << specific_heat(f, beta) << std::endl;
  }
}

int main(int argc, char **argv) {
  using namespace boost::multiprecision;
  options opt(argc, argv);
  if (!opt.valid) {
    std::cerr << "Usage: " << argv[0] << " [-p prec] T\n"
              << "       " << argv[0] << " [-p prec] J T\n"
              << "       " << argv[0] << " [-p prec] Ja Jb Jc T\n"
              << "       " << argv[0] << " [-p prec] Ja Jb Jc Tmin Tmax dT\n";
    return 127;
  }
  if (opt.prec <= std::numeric_limits<float>::digits10) {
    calc<float>(opt);
  } else if (opt.prec <= std::numeric_limits<double>::digits10) {
    calc<double>(opt);
  } else if (opt.prec <=
             std::numeric_limits<mp_wrapper<cpp_dec_float_50>>::digits10) {
    calc<mp_wrapper<cpp_dec_float_50>>(opt);
  } else if (opt.prec <=
             std::numeric_limits<mp_wrapper<cpp_dec_float_100>>::digits10) {
    calc<mp_wrapper<cpp_dec_float_100>>(opt);
  } else {
    std::cerr << "Error: Required precision is too high\n";
    return 127;
  }
}
