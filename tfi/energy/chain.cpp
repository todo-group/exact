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

// Ground state energy of the transverse Ising chain

#include <iomanip>
#include <iostream>
#include <string>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "options.hpp"
#include "chain.hpp"
#include "convert.hpp"

template<typename T>
void calc(const options& opt) {
  using namespace tfi::energy;
  using std::abs;
  typedef T real_t;
  real_t J = convert<real_t>(opt.J);
  real_t GammaMin, GammaMax, dGamma;
  if (opt.GammaMin == "GammaC") {
    GammaMin = GammaMax = dGamma = abs(J);
  } else {
    GammaMin = convert<real_t>(opt.GammaMin);
    GammaMax = convert<real_t>(opt.GammaMax);
    dGamma = convert<real_t>(opt.dGamma);
  }
  // if (Tmin < 0 || Tmax < 0) throw(std::invalid_argument("Temperature should be positive"));
  // if (GammaMin > GammaMax) throw(std::invalid_argument("GammaMax should be larger than GammaMin"));
  // if (dGamma <= 0) throw(std::invalid_argument("dGamma should be positive"));
  std::cout << std::scientific << std::setprecision(std::numeric_limits<real_t>::digits10)
            << "# lattice: chain\n"
            << "# precision: " << std::numeric_limits<real_t>::digits10 << std::endl
            << "# L J Gamma E/L\n";
  for (auto Gamma = GammaMin; Gamma < GammaMax + 1e-4 * dGamma; Gamma += dGamma) {
    auto ene = chain::infinite(J, Gamma);
    std::cout << "inf " << J << ' ' << Gamma << ' ' << ene << std::endl;
  }
}

int main(int argc, char **argv) {
  using namespace boost::multiprecision;
  options opt(argc, argv);
  if (!opt.valid) return 127;
  if (opt.prec <= std::numeric_limits<float>::digits10) {
    calc<float>(opt);
  } else if (opt.prec <= std::numeric_limits<double>::digits10) {
    calc<double>(opt);
  } else if (opt.prec <= std::numeric_limits<cpp_dec_float_50>::digits10) {
    calc<cpp_dec_float_50>(opt);
  } else if (opt.prec <= std::numeric_limits<cpp_dec_float_100>::digits10) {
    calc<cpp_dec_float_100>(opt);
  } else {
    std::cerr << "Error: Required precision is too high\n"; return 127;
  }
}
