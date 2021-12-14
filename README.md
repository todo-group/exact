# exact

* Numerically exact calculation of critical temperature of two-dimensional Ising models (ising/tc)
* Numerically exact calculation of free energy, energy, and specific heat of two-dimensional (finite-size/inifinit-size) Ising models (ising/free_energy)
* Numerically exact calculation of density of state of two-dimensional finite-size Ising model (ising/dos)
* Numerically exact calculation of free energy, gap, and entropy of quantum antiferromagnetic Heisenberg chain (afh/free_energy)
* Numerically exact calculation of ground-state energy of transverse-field Ising model (tfi/energy)

Multi-precision calculation is supported in ising/tc, ising/free_energy, tfi/energy.

## Prerequisites

* C++14 compiler
* CMake https://cmake.org
* Boost https://www.boost.org

## Compile & test

```
mkdir build
cd build
cmake ..
make
ctest
```

## Gallery

* gallery/ising-square-tc

  * calculation of the free energy, energy, and specific heat of the finite-size square Ising model at the critical temperature for L=2,4,8,...,65536 (with double precision and 50-digit precision)

## References

* Critical temperature of square lattice Ising model
  * L. Onsager, "Crystal Statistics. I. A Two-Dimensional Model with an Order-Disorder Transition ", [Phys. Rev. 65, 117 (1944)](https://doi.org/10.1103/PhysRev.65.117).
* Critical temperature of triangular lattice Ising model
  * J. Stephenson, "Ising‐Model Spin Correlations on the Triangular Lattice. IV. Anisotropic Ferromagnetic and Antiferromagnetic Lattices", [J. of Math. Phys. 11, 420 (1970)](https://doi.org/10.1063/1.1665155).

* Free energy of square lattice Ising model
  * B. Kastening, "Simplifying Kaufman’s solution of the two-dimensional Ising model", [Phys. Rev. E 64, 066106 (2001)](https://doi.org/10.1103/PhysRevE.64.066106).

* Free energy of triangular lattice Ising model
  * R. M. F. Houtappel, "Order-disorder in hexagonal lattices", [Physica 16, 425 (1950)](https://doi.org/10.1016/0031-8914(50)90130-3).

* Density of state of square lattice Ising model
  * P. Beale, "Exact Distribution of Energies in the Two-Dimensional Ising Model", [Phys. Rev. Lett. 76, 78-81 (1996)](https://doi.org/10.1103/PhysRevLett.76.78).

* Ground state energy of transverse field Ising model
  * G. B. Mbeng, A. Russomanno, and G. E. Santoro, "The quantum Ising chain for beginners", [arXiv:2009.09208](https://arxiv.org/abs/2009.09208).

## License

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
