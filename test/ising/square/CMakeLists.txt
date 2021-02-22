set(PROGS free_energy_finite counting_uniform counting_list counting_mag transfer_matrix transfer_matrix_uniform transfer_matrix_list finite)

foreach(name ${PROGS})
  set(target_name ising_square_${name})
  add_executable(${target_name} ${name}.cpp)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} ${ALPS_LIBRARIES})
  add_test(${name} ${name})
endforeach(name)
