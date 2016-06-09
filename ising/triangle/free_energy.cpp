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

// Calculating free energy density of triangle lattice Ising model

#include "free_energy.h"
#include <boost/lexical_cast.hpp>
#include <iostream>

int main(int argc, char **argv) {
  double t_min, t_max, t_step;
  int Nint;
  if (argc >=5) {
    t_min = boost::lexical_cast<double>(argv[1]);
    t_max = boost::lexical_cast<double>(argv[2]);
    t_step = boost::lexical_cast<double>(argv[3]);
    Nint = boost::lexical_cast<double>(argv[4]);
  } else {
    std::cin >> t_min >> t_max >> t_step >> Nint;
  }
  std::cout << "# Nint = " << Nint << std::endl;
  for (double t = t_min; t <= t_max; t += t_step) {
    double beta = 1 / t;
    std::cout << t << ' '
              << ising::triangle::free_energy_density(beta, 1, 1, 1, Nint) << std::endl;
  }
}
