free_energy.h: 無限系の三角格子イジング模型の自由エネルギーを計算する関数
  ising::triangle::free_energy_density(double beta, double J1, double J2, double J3, int Nint);

free_energy.cpp: サンプルプログラム
  使用法1) ./finite_energy Tmin Tmax Tstep Nint
  使用法2) ./finite_energy < parameter.txt

free_energy.op: ./finite_energy < free_energy.ip の出力結果

参考文献:
[1] R. M. F. Houtappel, Physica 16, 425 (1950).
