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

// Calculating free energy density of Ising model by exact counting

#include <stdexcept>
#include <tuple>
#include <vector>
#include <standards/exp_number.hpp>

#ifndef ISING_COUNTING_HPP
#define ISING_COUNTING_HPP

namespace ising {

struct counting {
public:
  template<typename LATTICE>
  static double free_energy(double beta, LATTICE const& lat, std::vector<double> const& inter,
                            std::vector<double> const& field = std::vector<double>(0)) {
    if (beta <= 0)
      throw(std::invalid_argument("beta should be positive"));
    if (lat.num_sites() > 30)
      throw(std::invalid_argument("too large lattice"));
    if (inter.size() != lat.num_bonds())
      throw(std::invalid_argument("inconsitent table size of interaction"));
    if (field.size() > 0 && field.size() != lat.num_sites())
      throw(std::invalid_argument("inconsitent table size of external field"));
    unsigned long num_states = 1 << lat.num_sites();
    standards::exp_double sum = 0;
    for (unsigned long c = 0; c < num_states; ++c) {
      standards::exp_double weight = 1;
      for (unsigned int b = 0; b < lat.num_bonds(); ++b) {
        int ci = (c >> lat.source(b)) & 1;
        int cj = (c >> lat.target(b)) & 1;
        weight *= standards::exp_number<double>(beta * inter[b] * (1 - 2 * (ci ^ cj)));
      }
      if (field.size()) {
        for (unsigned int s = 0; s < lat.num_sites(); ++s) {
          int c0 = (c >> s) & 1;
          weight *= standards::exp_number<double>(beta * field[s] * (1 - 2 * c0));
        }
      }
      sum += weight;
    }
    return -log(sum) / beta;
  }

  template<typename LATTICE>
  static std::tuple<double, double, double, double>
  magnetization(double beta, LATTICE const& lat, std::vector<double> const& inter,
                 std::vector<double> const& field = std::vector<double>(0)) {
    if (beta <= 0)
      throw(std::invalid_argument("beta should be positive"));
    if (lat.num_sites() > 30)
      throw(std::invalid_argument("too large lattice"));
    if (inter.size() != lat.num_bonds())
      throw(std::invalid_argument("inconsitent table size of interaction"));
    if (field.size() > 0 && field.size() != lat.num_sites())
      throw(std::invalid_argument("inconsitent table size of external field"));
    unsigned long num_states = 1 << lat.num_sites();
    standards::exp_double sum = 0;
    standards::exp_double sum_m1 = 0;
    standards::exp_double sum_m2 = 0;
    standards::exp_double sum_m3 = 0;
    standards::exp_double sum_m4 = 0;
    for (unsigned long c = 0; c < num_states; ++c) {
      double m = 0;
      for (unsigned int s = 0; s < lat.num_sites(); ++s) {
        m += 1.0- 2 * ((c >> s) & 1);
      }
      standards::exp_double weight = 1;
      for (unsigned int b = 0; b < lat.num_bonds(); ++b) {
        int ci = (c >> lat.source(b)) & 1;
        int cj = (c >> lat.target(b)) & 1;
        weight *= standards::exp_number<double>(beta * inter[b] * (1 - 2 * (ci ^ cj)));
      }
      if (field.size()) {
        for (unsigned int s = 0; s < lat.num_sites(); ++s) {
          int c0 = (c >> s) & 1;
          weight *= standards::exp_number<double>(beta * field[s] * (1 - 2 * c0));
        }
      }
      sum += weight;
      sum_m1 += std::pow(m, 1.0) * weight;
      sum_m2 += std::pow(m, 2.0) * weight;
      sum_m3 += std::pow(m, 3.0) * weight;
      sum_m4 += std::pow(m, 4.0) * weight;
    }
    return std::make_tuple(sum_m1 / sum, sum_m2 / sum, sum_m3 / sum, sum_m4 / sum);
  }

  template<typename LATTICE>
  static double free_energy(double beta, LATTICE const& lat, double J, double H = 0.0) {
    std::vector<double> inter(lat.num_bonds(), J);
    std::vector<double> field((H != 0.0) ? lat.num_sites() : 0, H);
    return free_energy(beta, lat, inter, field);
  }

  template<typename LATTICE>
  static std::tuple<double, double, double, double>
  magnetization(double beta, LATTICE const& lat, double J, double H = 0.0) {
    std::vector<double> inter(lat.num_bonds(), J);
    std::vector<double> field((H != 0.0) ? lat.num_sites() : 0, H);
    return magnetization(beta, lat, inter, field);
  }

  template<typename LATTICE>
  static double free_energy_density(double beta, LATTICE const& lat,
                                    std::vector<double> const& inter,
                                    std::vector<double> const& field = std::vector<double>(0)) {
    return free_energy(beta, lat, inter, field) / lat.num_sites();
  }

  template<typename LATTICE>
  static double free_energy_density(double beta, LATTICE const& lat, double J, double H = 0.0) {
    return free_energy(beta, lat, J, H) / lat.num_sites();
  }

  template<typename LATTICE>
  static std::tuple<double, double, double, double>
  magnetization_density(double beta, LATTICE const& lat,
                        std::vector<double> const& inter,
                        std::vector<double> const& field = std::vector<double>(0)) {
    double m1, m2, m3, m4;
    std::tie(m1, m2, m3, m4) = magnetization(beta, lat, inter, field);
    return std::make_tuple(m1 / lat.num_sites(), m2 / std::pow(lat.num_sites(), 2.0),
                             m3 / std::pow(lat.num_sites(), 3.0),
                             m4 / std::pow(lat.num_sites(), 4.0));
  }

  template<typename LATTICE>
  static std::tuple<double, double, double, double>
  magnetization_density(double beta, LATTICE const& lat, double J, double H = 0.0) {
    double m1, m2, m3, m4;
    std::tie(m1, m2, m3, m4) = magnetization(beta, lat, J, H);
    return std::make_tuple(m1 / lat.num_sites(), m2 / std::pow(lat.num_sites(), 2.0),
                             m3 / std::pow(lat.num_sites(), 3.0),
                             m4 / std::pow(lat.num_sites(), 4.0));
  }
};

} // end namespace ising

#endif // ISING_COUNTING_HPP
