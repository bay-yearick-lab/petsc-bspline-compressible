add_executable(test_operators test_operators.f90)
target_link_libraries(test_operators PRIVATE petb)

# Manufactured-solution checks (ctest), each on 4 MPI ranks.
set(PETB_TESTS linearsolve_x linearsolve_y linearsolve_z
               nonlinear_sum_BBB divergence viscous)
foreach(t ${PETB_TESTS})
  add_test(NAME mms_${t}
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                   ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_operators> -test ${t})
  set_tests_properties(mms_${t} PROPERTIES TIMEOUT 600)
endforeach()
