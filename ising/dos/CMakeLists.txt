set(PROGS square_count square_finite)

foreach(name ${PROGS})
  set(target_name ising_dos_${name})
  add_executable(${target_name} ${name}.cpp)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
endforeach(name)
