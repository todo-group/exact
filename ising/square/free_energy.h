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

// Calculating free energy density of square lattice Ising model in the thermodynamic limit

// reference: B. Kastening, Phys. Rev. E 64, 066106 (2001), wrn:2011/02/10

#ifndef ISING_SQUARE_FREE_ENERGY_H
#define ISING_SQUARE_FREE_ENERGY_H

#include <iostream>
#include <cmath>

namespace {

inline double f(double theta, double Ka, double Kb) {
  double k = 1.0 / (std::sinh(2 * Ka) * std::sinh(2 * Kb));
  return std::log(std::cosh(2 * Ka) * std::cosh(2 * Kb) +
                  std::sqrt(1 + k * k - 2 * k * std::cos(2 * theta)) / k);
}

}

namespace ising {
namespace square {
  
inline double free_energy_density(double beta, double Ja, double Jb, int n = 100) {
  double Ka = beta * Ja;
  double Kb = beta * Jb;
  double g = 0.0;
  double dt = M_PI / n;
  for (int i = 0; i < n; ++i) {
    double t0 = dt * i;
    double t1 = dt * (i+1);
    g += f(t0, Ka, Kb) + 4 * f((t0 + t1) / 2, Ka, Kb) + f(t1, Ka, Kb);
  }
  g *= dt / 6;
  g /= 2 * M_PI;
  return - (std::log(2.0) / 2 + g) / beta;
}

} // end namespace square
} // end namespace ising

#endif // ISING_SQUARE_FREE_ENERGY_H
