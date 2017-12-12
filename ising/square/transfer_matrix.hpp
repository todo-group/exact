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
  static lse::exp_double product_D(double beta, std::vector<double> const& interaction_x, VEC& v) {
    int width = interaction_x.size();
    int dim = 1 << width;
    lse::exp_double normal = 1;
    std::vector<boost::array<double, 2> > weight(width);
    for (int i = 0; i < width; ++i) {
      if (interaction_x[i] >= 0) {
        // ferromagnetic
        normal *= lse::exp_value(beta * interaction_x[i]);
        weight[i][0] = 1;
        weight[i][1] = std::exp(-2 * beta * interaction_x[i]);
      } else {
        // antiferromagnetic
        normal *= lse::exp_value(-beta * interaction_x[i]);
        weight[i][0] = std::exp(2 * beta * interaction_x[i]);
        weight[i][1] = 1;
      }
    }
    for (int s = 0; s < dim; ++s) {
      double elem = v[s];
      for (int i = 0; i < width; ++i) {
        int j = (i+1) % width;
        elem *= weight[i][((s >> i) & 1) ^ ((s >> j) & 1)];
      }
      v[s] = elem;
    }
    return normal;
  }

  template<typename VEC>
  static lse::exp_double product_U(double beta, std::vector<double> const& interaction_y, VEC& v) {
    int width = interaction_y.size();
    int dim = 1 << width;
    lse::exp_double normal = 1;
    boost::array<double, 2> weight;
    for (int i = 0; i < width; ++i) {
      if (interaction_y[i] >= 0) {
        // ferromagnetic
        normal *= lse::exp_value(beta * interaction_y[i]);
        weight[0] = 1;
        weight[1] = std::exp(-2 * beta * interaction_y[i]);
      } else {
        // antiferromagnetic
        normal *= lse::exp_value(-beta * interaction_y[i]);
        weight[0] = std::exp(2 * beta * interaction_y[i]);
        weight[1] = 1;
      }
      for (int s0 = 0; s0 < dim; ++s0) {
        if (((s0 >> i) & 1) == 0) {
          int s1 = s0 ^ (1 << i);
          double v0 = v[s0];
          double v1 = v[s1];
          v[s0] = weight[0] * v0 + weight[1] * v1;
          v[s1] = weight[0] * v1 + weight[1] * v0;
        }
      }
    }
    return normal;
  }
};

} // end namespace square
} // end namespace ising

#endif // ISING_SQUARE_TRANSFER_MATRIX_HPP
