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

// Critical temperature of hexagonal lattice Ising model

#include <iomanip>
#include <iostream>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "ising/mp_wrapper.hpp"
#include "options.hpp"
#include "hexagonal.hpp"

template<class T>
void calc(const std::string& Ja_in, const std::string& Jb_in, const std::string& Jc_in) {
  typedef T real_t;
  real_t Ja = convert<real_t>(Ja_in);
  real_t Jb = convert<real_t>(Jb_in);
  real_t Jc = convert<real_t>(Jc_in);
  auto tc = ising::tc::hexagonal(Ja, Jb, Jc);
  std::cout << std::scientific << std::setprecision(std::numeric_limits<real_t>::digits10)
            << "# lattice: hexagonal\n"
            << "# precision: " << std::numeric_limits<real_t>::digits10 << std::endl
            << "# Ja Jb Jc Tc 1/Tc" << std::endl
            << Ja << ' ' << Jb << ' ' << Jc << ' ' << tc << ' ' << (1 / tc) << std::endl;
}

int main(int argc, char **argv) {
  namespace mp = boost::multiprecision;
  options3 opt(argc, argv);
  if (!opt.valid) return 127;
  if (opt.prec <= std::numeric_limits<float>::digits10) {
    calc<float>(opt.Ja, opt.Jb, opt.Jc);
  } else if (opt.prec <= std::numeric_limits<double>::digits10) {
    calc<double>(opt.Ja, opt.Jb, opt.Jc);
  } else if (opt.prec <= std::numeric_limits<mp_wrapper<mp::cpp_dec_float_50>>::digits10) {
    calc<mp_wrapper<mp::cpp_dec_float_50>>(opt.Ja, opt.Jb, opt.Jc);
  } else if (opt.prec <= std::numeric_limits<mp_wrapper<mp::cpp_dec_float_100>>::digits10) {
    calc<mp_wrapper<mp::cpp_dec_float_100>>(opt.Ja, opt.Jb, opt.Jc);
  } else {
    std::cerr << "Error: Required precision is too high\n"; return 127;
  }
}
