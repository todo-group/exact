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

#include "free_energy_finite.h"
#include <boost/lexical_cast.hpp>

int main(int argc, char **argv) {
  int L; // system size
  double t_min, t_max, t_step;
  if (argc >=5) {
    L = boost::lexical_cast<int>(argv[1]);
    t_min = boost::lexical_cast<double>(argv[2]);
    t_max = boost::lexical_cast<double>(argv[3]);
    t_step = boost::lexical_cast<double>(argv[4]);
  } else {
    std::cin >> L >> t_min >> t_max >> t_step;
  }
  std::cout << "# L = " << L << std::endl;
  for (double t = t_min; t <= t_max; t += t_step) {
    double beta = 1 / t;
    std::cout << t << ' '
              << ising::square::free_energy_density(beta, -1, -1, L, L) << std::endl;
  }
}
