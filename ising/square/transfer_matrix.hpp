/*****************************************************************************
*
* Copyright (C) 2015-2017 by Synge Todo <wistaria@phy.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

// Calculating free energy density of square lattice Ising model by the transfer matrix method

#include <lattice/square.hpp>
#include <lse/exp_number.hpp>
#include <boost/array.hpp>
#include <cmath>
#include <vector>

#ifndef ISING_SQUARE_TRANSFER_MATRIX_HPP
#define ISING_SQUARE_TRANSFER_MATRIX_HPP

namespace ising {
namespace square {

struct transfer_matrix {
public:
  template<typename VEC>
  static lse::exp_double product_D(double beta, std::vector<double> const& inter_x,
                                   std::vector<double> const& field, VEC& v) {
    int width = inter_x.size();
    int dim = 1 << width;
    lse::exp_double normal = 1;
    std::vector<boost::array<double, 2> > weight(width);
    for (int b = 0; b < width; ++b) {
      double offset = std::abs(beta * inter_x[b]);
      normal *= lse::exp_value(offset);
      weight[b][0] = std::exp(beta * inter_x[b] - offset);
      weight[b][1] = std::exp(-beta * inter_x[b] - offset);
    }
    for (int c = 0; c < dim; ++c) {
      double elem = v[c];
      for (int b = 0; b < width; ++b) {
        int s0 = b;
        int s1 = (b+1) % width;
        elem *= weight[b][((c >> s0) & 1) ^ ((c >> s1) & 1)];
      }
      v[c] = elem;
    }
    if (field.size()) {
      for (int s = 0; s < width; ++s) {
        double offset = std::abs(beta * field[s]);
        normal *= lse::exp_value(offset);
        weight[s][0] = std::exp(beta * field[s] - offset);
        weight[s][1] = std::exp(-beta * field[s] - offset);
      }
      for (int c = 0; c < dim; ++c) {
        double elem = v[c];
        for (int s = 0; s < width; ++s) elem *= weight[s][((c >> s) & 1)];
        v[c] = elem;
      }
    }
    return normal;
  }

  template<typename VEC>
  static lse::exp_double product_D(double beta, std::vector<double> const& inter_x, VEC& v) {
    std::vector<double> field(0);
    return product_D(beta, inter_x, field, v);
  }

  template<typename VEC>
  static lse::exp_double product_U(double beta, std::vector<double> const& inter_y, VEC& v) {
    int width = inter_y.size();
    int dim = 1 << width;
    lse::exp_double normal = 1;
    boost::array<double, 2> weight;
    for (int s = 0; s < width; ++s) {
      double offset = std::abs(beta * inter_y[s]);
      normal *= lse::exp_value(offset);
      weight[0] = std::exp(beta * inter_y[s] - offset);
      weight[1] = std::exp(-beta * inter_y[s] - offset);
      for (int c0 = 0; c0 < dim; ++c0) {
        if (((c0 >> s) & 1) == 0) {
          int c1 = c0 ^ (1 << s);
          double v0 = v[c0];
          double v1 = v[c1];
          v[c0] = weight[0] * v0 + weight[1] * v1;
          v[c1] = weight[0] * v1 + weight[1] * v0;
        }
      }
    }
    return normal;
  }

  static double free_energy(double beta, lattice::square const& lat,
                            std::vector<double> const& inter,
                            std::vector<double> const& field = std::vector<double>(0)) {
    int Lx = lat.get_length_x();
    int Ly = lat.get_length_y();
    std::vector<double> inter_x(Lx), inter_y(Lx), field_x(field.size() > 0 ? Lx : 0);
    int dim = 1 << Lx;
    std::vector<double> v(dim);
    lse::exp_double sum = 0;
    for (int i = 0; i < dim; ++i) {
      lse::exp_double weight = 1;
      for (int j = 0; j < dim; ++j) v[j] = 0;
      v[i] = 1;
      for (int y = 0; y < Ly; ++y) {
        for (int x = 0; x < Lx; ++x) {
          inter_x[x] = inter[2 * (Lx * y + x)];
          inter_y[x] = inter[2 * (Lx * y + x) + 1];
        }
        if (field_x.size()) {
          for (int x = 0; x < Lx; ++x) {
            field_x[x] = field[Lx * y + x];
          }
        }
        weight *= ising::square::transfer_matrix::product_D(beta, inter_x, field_x, v);
        weight *= ising::square::transfer_matrix::product_U(beta, inter_y, v);
      }
      weight *= v[i];
      sum += weight;
    }
    return -log(sum) / beta;
  }

  static double free_energy(double beta, lattice::square const& lat, double J, double H = 0.0) {
    std::vector<double> inter(lat.num_bonds(), J);
    std::vector<double> field((H != 0.0) ? lat.num_sites() : 0, H);
    return free_energy(beta, lat, inter, field);
  }

  static double free_energy_density(double beta, lattice::square const& lat,
                                    std::vector<double> const& inter,
                                    std::vector<double> const& field = std::vector<double>(0)) {
    return free_energy(beta, lat, inter, field) / lat.num_sites();
  }

  static double free_energy_density(double beta, lattice::square const& lat, double J,
                                    double H = 0.0) {
    return free_energy(beta, lat, J, H) / lat.num_sites();
  }
};

} // end namespace square
} // end namespace ising

#endif // ISING_SQUARE_TRANSFER_MATRIX_HPP
