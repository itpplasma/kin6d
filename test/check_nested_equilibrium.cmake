file(MAKE_DIRECTORY "${WORK}/output")
execute_process(COMMAND "${APP}" --config "${CONFIG}" WORKING_DIRECTORY "${WORK}"
    RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE errors)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "nested_equilibrium failed: ${output}${errors}")
endif()
separate_arguments(expected UNIX_COMMAND "${EXPECTED}")
execute_process(COMMAND "${CHECK}" "${WORK}/output/nested_equilibrium.nc" ${expected}
    RESULT_VARIABLE status OUTPUT_VARIABLE check_output ERROR_VARIABLE errors)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "result check failed: ${check_output}${errors}")
endif()
message(STATUS "${check_output}")
