# test/ising/square

Example/test programs for exact free energy calculation of
square-lattice Ising model

## uniform (J=1) ferromagnetic Ising model without magnetic field

### free\_energy: free energy density in the thermodynamic limit

```
./free_energy t_min t_max t_step Nint
```

* t\_min, t\_max, t\_step: minimum/maximum/interval of temperature
* Nint: number of mesh in numerical integration
* calculation cost: Nint * (t\_max - t\_min) / t\_step
* memory cost: O(1)

### free\_energy\_finite: free energy density of LxL square lattice

```
./free_energy_finite L t_min t_max t_step
```

* L: linear size of lattice
* t\_min, t\_max, t\_step: minimum/maximum/interval of temperature
* calculation cost: L * (t\_max - t\_min) / t\_step
* memory cost: O(1)

## uniform Ising model in external field

### counting\_uniform: free energy density by exact counting

```
./counting_uniform L J H t_min t_max t_step
```
* L: linear size of lattice
* J: (uniform) coupling constant
* H: (uniform) external field
* t\_min, t\_max, t\_step: minimum/maximum/interval of temperature
* calculation cost: (L*L) * 2^(L*L) * (t\_max - t\_min) / t\_step
* memory cost: O(1)

### transfer\_matrix\_uniform: free energy density by transfer matrix method

```
./transfer_matrix_uniform L J H t_min t_max t_step
```
* L: linear size of lattice
* J: (uniform) coupling constant
* H: (uniform) external field
* t\_min, t\_max, t\_step: minimum/maximum/interval of temperature
* calculation cost: L * 2^L * (t\_max - t\_min) / t\_step
* memory cost: 2^L

## Ising model with random coupling and random field

### counting\_list: free energy density by exact counting

```
./counting_list < parameter_list
```

* format of parameter_list file

  ```
  Lx Ly
  J[0] J[1] ... J[2*Lx*Ly-1]
  H[0] H[1] ... H[Lx*Ly-1]
  t_min t_max t_step
  ```
* Lx, Ly: linear sizes of lattice in x- and y-directions
* J[b]: coupling constant of b-th bond
* H[s]: external field at s-th site
* t\_min, t\_max, t\_step: minimum/maximum/interval of temperature
* calculation cost: (L*L) * 2^(L*L) * (t\_max - t\_min) / t\_step
* memory cost: O(1)

### transfer\_matrix\_list: free energy density by transfer matrix method

```
./transfer_matrix_list < parameter_list
```

* format of parameter_list file

  ```
  Lx Ly
  J[0] J[1] ... J[2*Lx*Ly-1]
  H[0] H[1] ... H[Lx*Ly-1]
  t_min t_max t_step
  ```
* Lx, Ly: linear sizes of lattice in x- and y-directions
* J[b]: coupling constant of b-th bond
* H[s]: external field at s-th site
* t\_min, t\_max, t\_step: minimum/maximum/interval of temperature
* calculation cost: L * 2^L * (t_max - t_min) / t_step
* memory cost: 2^L
