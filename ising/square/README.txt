free_energy_finite.h: 有限系の2次元イジング模型の自由エネルギーを計算する関数
  ising::square::partition_function(double beta, double Ja, double Jb, int n, int m);
  ising::square::free_energy(double beta, double Ja, double Jb, int n, int m);
  ising::square::free_energy_density(double beta, double Ja, double Jb, int n, int m);

free_energy_finite.cpp: サンプルプログラム
  使用法1) ./finite_energy_finite Lx Ly Tmin Tmax Tstep
  使用法2) ./finite_energy_finite < parameter.txt

free_energy_fintie.op: ./finite_energy_finite < free_energy_finite.ip の出力結果

参考文献:
[1] B. Kastening, Phys. Rev. E 64, 066106 (2001).
