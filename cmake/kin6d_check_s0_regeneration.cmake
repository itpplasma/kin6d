file(MAKE_DIRECTORY "${WORK}/src/generated")
execute_process(COMMAND "${GENERATOR}" WORKING_DIRECTORY "${WORK}"
    RESULT_VARIABLE result OUTPUT_VARIABLE evidence ERROR_VARIABLE errors)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "gen/gen_s0.f90: ${evidence}${errors}")
endif()
execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files "${EXPECTED}"
    "${WORK}/src/generated/kin6d_s0_leaves.f90" RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "src/generated/kin6d_s0_leaves.f90: regeneration differs")
endif()
message(STATUS "${evidence}Regeneration matches the checked-in output")
