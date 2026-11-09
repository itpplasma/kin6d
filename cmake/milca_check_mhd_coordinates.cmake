file(MAKE_DIRECTORY "${WORK}")
execute_process(COMMAND "${GENERATOR}" "${WORK}" "${WORK}/identities.json" "${REVISION}"
    RESULT_VARIABLE result OUTPUT_VARIABLE evidence ERROR_VARIABLE errors)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "Coordinate identities: ${evidence}${errors}")
endif()
execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files "${SOURCE}"
    "${WORK}/ideal_mhd_energy_generated.f90" RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "Shared ideal-MHD coordinate regeneration differs")
endif()
message(STATUS "${evidence}Shared ideal-MHD leaves reproduced")
