set(PROGS square square_finite square_finite_tc square_count square_count_tc triangular)

foreach(name ${PROGS})
  set(target_name ising_free_energy_${name})
  add_executable(${target_name} ${name}.cpp)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
endforeach(name)
