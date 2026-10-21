# Public synthetic regression; historical CTest name retained for consumers.
# The fixture README records its independent analytic oracle and recipe.
# The three estimator fields now reference the RT equilibration introduced
# by 5dc6103, rather than its functional-recovery predecessor. The physical
# state, profiles and field references remain the M0 values. Check the method
# and error scope as well as the numbers: this nonlinear case has only a
# sampled frozen-source residual estimate, not a continuous certificate.
# Exact-error coverage lives in grad_shafranov_equilibrated_flux_majorant.
file(REMOVE_RECURSE "${WORK}")
file(MAKE_DIRECTORY "${WORK}")
foreach(name boundary.dat source.dat input.nml points.dat)
    file(COPY "${FIXTURES}/${name}" DESTINATION "${WORK}")
endforeach()
execute_process(COMMAND "${APP}" --config input.nml WORKING_DIRECTORY "${WORK}"
    OUTPUT_FILE "${WORK}/native.log" ERROR_FILE "${WORK}/native.err"
    RESULT_VARIABLE status)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "K6-mini gs_phase1 failed (${status}); see ${WORK}/native.log")
endif()
file(READ "${WORK}/summary.dat" summary)
if(NOT summary MATCHES "majorant_method equilibrated-rt" OR
   NOT summary MATCHES "majorant_quadrature_status approximate" OR
   NOT summary MATCHES "majorant_scope frozen-source-energy-error" OR
   NOT summary MATCHES "majorant_status 4")
    message(FATAL_ERROR "K6-mini estimator method or error scope changed")
endif()
execute_process(COMMAND "${APP}" --evaluate mesh.dat 6 0 profile-law.dat
    WORKING_DIRECTORY "${WORK}" INPUT_FILE "${WORK}/points.dat"
    OUTPUT_FILE "${WORK}/field.dat" RESULT_VARIABLE status)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "K6-mini field query failed (${status})")
endif()
execute_process(COMMAND "${CHECK}" compare "${WORK}" "${FIXTURES}"
    RESULT_VARIABLE status)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "K6-mini outputs differ from the fixture")
endif()
