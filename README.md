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
