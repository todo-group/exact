set(PROGS free_energy)

foreach(name ${PROGS})
  set(target_name ising_triangle_${name})
  add_executable(${target_name} ${name}.cpp)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} ${ALPS_LIBRARIES})
  add_test(${name} ${name})
endforeach(name)
