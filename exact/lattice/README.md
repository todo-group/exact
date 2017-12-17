# Numbering convention in lattice classes

## lattice::square

### site numbering

```
+-> x      0        -- 1    -- ... -- Lx-1
|          Lx       -- Lx+1 -- ... -- 2Lx-1
v          ...
y          Lx(Ly-1) --    ...     --- LxLy-1
```

### neighboring sites

```
+-> x                   neighbor(s,3)
|                            |
v        neighbor(s,2) -- site s -- neighbor(s,0)
y                            |
                        neighbor(s,1)
```

### bond numbering

Up (left) site of the bond is called *source*, and
down (right) site is called *target*

```
+-> x                   bond 2(s-Lx)+1
|                            |
v          bond 2(s-1) -- site s -- bond 2s
y                            |
                        bond 2s+1
```

### example (Lx = 4, Ly = 3)

* num_sites = 12
* num_bonds = 24

```
          [ 8]    [ 9]    [10]    [11]
            |       |       |       |
           17      19      21       23
            |       |       |       |
  [ 3]--6-[ 0]--0-[ 1]--2-[ 2]--4-[ 3]--6-[ 0]
            |       |       |       |
            1       3       5       7
            |       |       |       |
  [ 7]-14-[ 4]--8-[ 5]-10-[ 6]-12-[ 7]-14-[ 4]
            |       |       |       |
            9      11      13      15
            |       |       |       |
  [11]-22-[ 8]-16-[ 9]-18-[10]-20-[11]-22-[ 8]
            |       |       |       |
           17      19      21      23
            |       |       |       |
          [ 0]    [ 1]    [ 2]    [ 3]
```