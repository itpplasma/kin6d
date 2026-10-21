file(MAKE_DIRECTORY "${WORK}")
execute_process(COMMAND "${APP}" --config "${CONFIG}" WORKING_DIRECTORY "${WORK}"
    RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE errors)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "nested_oracle_export failed: ${output}${errors}")
endif()
execute_process(COMMAND "${CHECK}" "${WORK}/output/input.nested_oracle" "${FAMILY}"
    WORKING_DIRECTORY "${WORK}" RESULT_VARIABLE status
    OUTPUT_VARIABLE check_output ERROR_VARIABLE errors)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "boundary consumer failed: ${check_output}${errors}")
endif()
message(STATUS "${output}${check_output}")
if(FAMILY EQUAL 1)
    execute_process(COMMAND "${CHECK}" "${WORK}/output/input.nested_equilibrium"
        "${FAMILY}" matched WORKING_DIRECTORY "${WORK}" RESULT_VARIABLE status
        OUTPUT_VARIABLE check_output ERROR_VARIABLE errors)
    if(NOT status EQUAL 0)
        message(FATAL_ERROR "matched input failed: ${check_output}${errors}")
    endif()
    message(STATUS "${check_output}")
endif()
