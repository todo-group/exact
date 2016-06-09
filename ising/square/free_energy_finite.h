/*****************************************************************************
*
* Copyright (C) 2011-2016 by Synge Todo <wistaria@comp-phys.org>
*
* This software is published under the ALPS Application License; you
* can use, redistribute it and/or modify it under the terms of the
* license, either version 1 or (at your option) any later version.
* 
* You should have received a copy of the ALPS Application License
* along with this software; see the file LICENSE. If not, the license
* is also available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

// Calculating free energy density of square lattice Ising model

// reference: B. Kastening, Phys. Rev. E 64, 066106 (2001), wrn:2011/02/10

#ifndef ISING_SQUARE_FREE_ENERGY_FINITE_H
#define ISING_SQUARE_FREE_ENERGY_FINITE_H

#include <alps/parapack/exp_number.h>
#include <vector>

namespace {
  
inline alps::exp_double cosh_value(double x) {
  return (alps::exp_value(x) + alps::exp_value(-x)) / 2;
}

inline alps::exp_double sinh_value(double x) {
  return (alps::exp_value(x) - alps::exp_value(-x)) / 2;
}

}

namespace ising {
namespace square {
  
inline alps::exp_double partition_function(double beta, double Ja, double Jb, int n, int m) {
  double a = - beta * Ja;
  double b = - beta * Jb;
  std::vector<double> gamma(2 * n);
  for (int k = 0; k < 2 * n; ++k) {
    alps::exp_double cosh_g =
      (cosh_value(2*a) * cosh_value(2*b) - cos(M_PI*k/n) * sinh_value(2*b)) / sinh_value(2*a);
    gamma[k] = log(cosh_g + sqrt(cosh_g * cosh_g - 1));
  }
  if (sinh_value(2*a) * sinh_value(2*b) > 1) gamma[0] = -gamma[0];
  alps::exp_double p0(1), p1(1), p2(1), p3(1);
  for (int k = 1; k <= n; ++k) {
    p0 *= 2 * cosh_value(m * gamma[2*k-1] / 2);
    p1 *= 2 * sinh_value(m * gamma[2*k-1] / 2);
    p2 *= 2 * cosh_value(m * gamma[2*k-2] / 2);
    p3 *= 2 * sinh_value(m * gamma[2*k-2] / 2);
  }
  alps::exp_double z = 0.5 * pow(2 * sinh_value(2*a), m*n/2) * (p0 + p1 + p2 - p3);
  return z;
}

inline double free_energy(double beta, double Ja, double Jb, int n, int m) {
  return -log(partition_function(beta, Ja, Jb, n, m)) / beta;
}

inline double free_energy_density(double beta, double Ja, double Jb, int n, int m) {
  return free_energy(beta, Ja, Jb, n, m) / (n*m);
}

} // end namespace square
} // end namespace ising

#endif // ISING_SQUARE_FREE_ENERGY_FINITE_H
