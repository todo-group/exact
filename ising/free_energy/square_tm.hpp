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

// Free energy, energy, specific heat, magnetization, magnetization^2 of square lattice Ising model
// Transfer matrix method

#pragma once

#include <cmath>
#include <stdexcept>
#include <boost/array.hpp>
#include <standards/exp_number.hpp>

namespace ising {
namespace free_energy {
namespace square {

template<typename T>
struct transfer_matrix {
public:
  typedef std::size_t uint_t;
  typedef T real_t;
  typedef standards::exp_number<real_t> exp_number;

  struct result_t {
    typedef real_t root_type;
    void set(uint_t i, real_t value) {
      derivatives_[i][0] = value;
    }
    void set(uint_t i, uint_t j, real_t value) {
      derivatives_[i][j] = value;
    }
    real_t derivative(uint_t i) {
      return derivatives_[i][0];
    }
    real_t derivative(uint_t i, uint_t j) {
      return derivatives_[i][j];
    }
    real_t derivatives_[3][3];
  };

  static void fill_D(uint_t Ly, real_t Jy, real_t beta, real_t h, std::vector<real_t>& diag, std::vector<real_t>& diag_10, std::vector<real_t>& diag_20, std::vector<real_t>& diag_01, std::vector<real_t>& diag_02) {
    uint_t dim = 1 << Ly;
    diag.resize(dim);
    diag_10.resize(dim);
    diag_20.resize(dim);
    diag_01.resize(dim);
    diag_02.resize(dim);
    for (uint_t c = 0; c < dim; ++c) {
      diag[c] = 1;
      diag_10[c] = 1;
      diag_20[c] = 1;
      diag_01[c] = 1;
      diag_02[c] = 1;
    }
    boost::array<real_t, 2> weight;
    weight[0] = std::exp(beta * Jy);
    weight[1] = std::exp(-beta * Jy);
    for (uint_t c = 0; c < dim; ++c) {
      real_t elem = diag[c];
      for (uint_t b = 0; b < Ly; ++b) {
        uint_t s0 = b;
        uint_t s1 = (b+1) % Ly;
        elem *= weight[((c >> s0) & 1) ^ ((c >> s1) & 1)];
      }
      diag[c] = elem;
    }
    weight[0] = std::exp(beta * h);
    weight[1] = std::exp(-beta * h);
    for (uint_t c = 0; c < dim; ++c) {
      real_t elem = diag[c];
      for (uint_t s = 0; s < Ly; ++s) {
        elem *= weight[((c >> s) & 1)];
      }
      diag[c] = elem;
    }
    for (uint_t c = 0; c < dim; ++c) {
      real_t factor = 0;
      for (uint_t b = 0; b < Ly; ++b) {
        uint_t s0 = b;
        uint_t s1 = (b+1) % Ly;
        factor += Jy * (-2.0 * ((c >> s0) & 1) + 1) * (-2.0 * ((c >> s1) & 1) + 1);
      }
      for (uint_t s = 0; s < Ly; ++s) {
        factor += h * (-2.0 * ((c >> s) & 1) + 1);
      }
      diag_10[c] = diag[c] * factor;
      diag_20[c] = diag[c] * factor * factor;
    }
    for (uint_t c = 0; c < dim; ++c) {
      real_t factor = 0;
      for (uint_t s = 0; s < Ly; ++s) {
        factor += beta * (-2.0 * ((c >> s) & 1) + 1);
      }
      diag_01[c] = diag[c] * factor;
      diag_02[c] = diag[c] * factor * factor;
    }
  }

  static exp_number product_D(std::vector<real_t> const& diag, std::vector<real_t> const& diag_10, std::vector<real_t> const& diag_20, std::vector<real_t> const& diag_01, std::vector<real_t> const& diag_02, exp_number factor, std::vector<real_t>& v, std::vector<real_t>& v_10, std::vector<real_t>& v_20, std::vector<real_t>& v_01, std::vector<real_t>& v_02) {
    uint_t dim = diag.size();
    real_t norm2 = 0;
    for (uint_t c = 0; c < dim; ++c) {
      v_20[c] = diag_20[c] * v[c] + 2 * diag_10[c] * v_10[c] + diag[c] * v_20[c];
      v_10[c] = diag_10[c] * v[c] + diag[c] * v_10[c];
      v_02[c] = diag_02[c] * v[c] + 2 * diag_01[c] * v_01[c] + diag[c] * v_02[c];
      v_01[c] = diag_01[c] * v[c] + diag[c] * v_01[c];
      v[c] = diag[c] * v[c];
      norm2 += v[c] * v[c];
    }
    auto norm_inv = 1 / std::sqrt(norm2);
    for (uint_t c = 0; c < dim; ++c) {
      v_20[c] *= norm_inv;
      v_10[c] *= norm_inv;
      v_02[c] *= norm_inv;
      v_01[c] *= norm_inv;
      v[c] *= norm_inv;
    }
    return factor / norm_inv;
  }

  static exp_number product_U(uint_t Ly, real_t Jx, real_t beta, exp_number factor, std::vector<real_t>& v, std::vector<real_t>& v_10, std::vector<real_t>& v_20, std::vector<real_t>& v_01, std::vector<real_t>& v_02, std::vector<real_t>& vt) {
    uint_t dim = 1 << Ly;
    boost::array<double, 2> weight;
    weight[0] = std::exp(beta * Jx);
    weight[1] = std::exp(-beta * Jx);
    for (uint_t s = 0; s < Ly; ++s) {
      for (uint_t c0 = 0; c0 < dim; ++c0) {
        uint_t c1 = c0 ^ (1 << s);
        vt[c0] = Jx * Jx * (weight[0] * v[c0] + weight[1] * v[c1]) + 2 * Jx * (weight[0] * v_10[c0] - weight[1] * v_10[c1]) + (weight[0] * v_20[c0] + weight[1] * v_20[c1]);
      }
      std::swap(v_20, vt);
      for (uint_t c0 = 0; c0 < dim; ++c0) {
        uint_t c1 = c0 ^ (1 << s);
        vt[c0] = Jx * (weight[0] * v[c0] - weight[1] * v[c1]) + (weight[0] * v_10[c0] + weight[1] * v_10[c1]);
      }
      std::swap(v_10, vt);
      for (uint_t c0 = 0; c0 < dim; ++c0) {
        uint_t c1 = c0 ^ (1 << s);
        vt[c0] = (weight[0] * v_02[c0] + weight[1] * v_02[c1]);
      }
      std::swap(v_02, vt);
      for (uint_t c0 = 0; c0 < dim; ++c0) {
        uint_t c1 = c0 ^ (1 << s);
        vt[c0] = (weight[0] * v_01[c0] + weight[1] * v_01[c1]);
      }
      std::swap(v_01, vt);
      for (uint_t c0 = 0; c0 < dim; ++c0) {
        uint_t c1 = c0 ^ (1 << s);
        vt[c0] = weight[0] * v[c0] + weight[1] * v[c1];
      }
      std::swap(v, vt);
    }
    real_t norm2 = 0;
    for (uint_t c = 0; c < dim; ++c) {
      norm2 += v[c] * v[c];
    }
    real_t norm_inv = 1 / std::sqrt(norm2);
    for (uint_t c = 0; c < dim; ++c) {
      v_20[c] *= norm_inv;
      v_10[c] *= norm_inv;
      v_02[c] *= norm_inv;
      v_01[c] *= norm_inv;
      v[c] *= norm_inv;
    }
    return factor / norm_inv;
  }
  
  static result_t calc(uint_t Lx, uint_t Ly, real_t Jx, real_t Jy, real_t beta, real_t h) {
    if (Lx == 0 || Ly == 0) throw(std::invalid_argument("Lx and Ly should be positive"));
    if (beta <= 0) throw(std::invalid_argument("beta should be positive"));
    uint_t dim = 1 << Ly;

    std::vector<real_t> diag(dim), diag_10(dim), diag_20(dim), diag_01(dim), diag_02(dim);
    fill_D(Ly, Jy, beta, h, diag, diag_10, diag_20, diag_01, diag_02);

    result_t res;
    exp_number sum = 0;
    exp_number sum_10 = 0;
    exp_number sum_20 = 0;
    exp_number sum_01 = 0;
    exp_number sum_02 = 0;
    std::vector<real_t> v(dim), v_10(dim), v_20(dim), v_01(dim), v_02(dim), vt(dim);
    for (uint_t i = 0; i < dim; ++i) {
      exp_number factor = 1;
      for (uint_t j = 0; j < dim; ++j) {
        v[j] = 0;
        v_10[j] = 0;
        v_20[j] = 0;
        v_01[j] = 0;
        v_02[j] = 0;
      }
      v[i] = 1;
      for (uint_t x = 0; x < Lx; ++x) {
        factor = transfer_matrix<real_t>::product_D(diag, diag_10, diag_20, diag_01, diag_02, factor, v, v_10, v_20, v_01, v_02);
        factor = transfer_matrix<real_t>::product_U(Ly, Jx, beta, factor, v, v_10, v_20, v_01, v_02, vt);
      }
      sum += factor * v[i];
      sum_10 += factor * v_10[i];
      sum_20 += factor * v_20[i];
      sum_01 += factor * v_01[i];
      sum_02 += factor * v_02[i];
    }
    res.set(0, 0, -log(sum) / (Lx * Ly * beta));
    res.set(1, 0, ((log(sum) / beta) - (sum_10 / sum)) / (Lx * Ly * beta));
    res.set(2, 0, (-(2 * log(sum) / beta / beta) + 2 * (sum_10 / sum / beta) - ((sum_20 * sum - sum_10 * sum_10) / sum / sum)) / (Lx * Ly * beta));
    res.set(0, 1, -sum_01 / sum / (Lx * Ly * beta));
    res.set(0, 2, -sum_02 / sum / (Lx * Ly * beta));
    return res;
  }
};

} // end namespace square
} // end namespace free_energy
} // end namespace ising
