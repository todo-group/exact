set(PROGS square triangular hexagonal)

foreach(name ${PROGS})
  set(target_name test_ising_tc_${name})
  add_executable(${target_name} ${name}.cpp)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} gtest_main)
  add_test(${target_name} ${name})
endforeach(name)
