# Run GENERATOR in WORK and compare src/generated/LEAF with EXPECTED.
file(MAKE_DIRECTORY "${WORK}/src/generated")
execute_process(COMMAND "${GENERATOR}" WORKING_DIRECTORY "${WORK}"
    RESULT_VARIABLE result OUTPUT_VARIABLE evidence ERROR_VARIABLE errors)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "${GENERATOR}: ${evidence}${errors}")
endif()
execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files "${EXPECTED}"
    "${WORK}/src/generated/${LEAF}" RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "src/generated/${LEAF}: regeneration differs")
endif()
message(STATUS "${evidence}Regeneration matches the checked-in output")
