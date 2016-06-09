set(PROGS free_energy free_energy_finite)

foreach(name ${PROGS})
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} ${ALPS_LIBRARIES})
  add_alps_test(${name})
endforeach(name)
