free_energy.h: 無限系の2次元イジング模型の自由エネルギーを計算する関数
  ising::square::free_energy_density(double beta, double Jx, double Jy, int Nint);

free_energy_finite.h: 有限系の2次元イジング模型の自由エネルギーを計算する関数
  ising::square::partition_function(double beta, double Jx, double Jy, int Lx, int Ly);
  ising::square::free_energy(double beta, double Jx, double Jy, int Lx, int Ly);
  ising::square::free_energy_density(double beta, double Jx, double Jy, int Lx, int Ly);

free_energy.cpp: サンプルプログラム
  使用法1) ./finite_energy Tmin Tmax Tstep Nint
  使用法2) ./finite_energy < parameter.txt

free_energy_finite.cpp: サンプルプログラム
  使用法1) ./finite_energy_finite L Tmin Tmax Tstep
  使用法2) ./finite_energy_finite < parameter.txt

free_energy.op: ./finite_energy < free_energy.ip の出力結果

free_energy_fintie.op: ./finite_energy_finite < free_energy_finite.ip の出力結果

参考文献:
[1] https://en.wikipedia.org/wiki/Square-lattice_Ising_model
[2] B. Kastening, Phys. Rev. E 64, 066106 (2001).
