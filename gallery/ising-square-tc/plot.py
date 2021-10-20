import math
import matplotlib.pyplot as plt

filename = "result-p15.dat"

with open(filename, 'r') as f:
    for line in f:
        data = line.split()
        if (data[0] == "inf"):
            free_energy_inf = float(data[6])
            energy_inf = float(data[7])

L = []
free_energy = []
energy = []
with open(filename, 'r') as f:
    for line in f:
        data = line.split()
        if (data[0] != '#' and data[0] != "inf"):
            L.append(int(data[0]))
            free_energy.append(math.fabs(float(data[6])-free_energy_inf))
            energy.append(math.fabs(float(data[7])-energy_inf))
            
plt.plot(L, free_energy, marker = 'o', label = 'finite-size error of free energy density')
plt.plot(L, energy, marker = 'v', label = 'finite-size error of energy density')
plt.xlabel('L')
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.legend()
plt.savefig('plot.pdf')
plt.show()

print(L)
print(free_energy)
print(energy)
