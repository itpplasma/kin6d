# Run GENERATOR in WORK and compare GENERATED_DIR/LEAF with EXPECTED.
if(NOT DEFINED GENERATED_DIR)
    set(GENERATED_DIR src/generated)
endif()
file(MAKE_DIRECTORY "${WORK}/${GENERATED_DIR}")
execute_process(COMMAND "${GENERATOR}" WORKING_DIRECTORY "${WORK}"
    RESULT_VARIABLE result OUTPUT_VARIABLE evidence ERROR_VARIABLE errors)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "${GENERATOR}: ${evidence}${errors}")
endif()
execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files "${EXPECTED}"
    "${WORK}/${GENERATED_DIR}/${LEAF}" RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "${GENERATED_DIR}/${LEAF}: regeneration differs")
endif()
message(STATUS "${evidence}Regeneration matches the checked-in output")
