file(MAKE_DIRECTORY "${WORK}/src/generated")
execute_process(COMMAND "${GENERATOR}" WORKING_DIRECTORY "${WORK}"
    RESULT_VARIABLE result OUTPUT_VARIABLE evidence ERROR_VARIABLE errors)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "gen/gen_vrmhd_leaves.f90: ${evidence}${errors}")
endif()
execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files "${EXPECTED}"
    "${WORK}/src/generated/milca_vrmhd_leaves.f90" RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "milca_vrmhd_leaves.f90: regeneration differs")
endif()

get_filename_component(expected_directory "${EXPECTED}" DIRECTORY)
execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${expected_directory}/milca_vrmhd_dispatch.f90"
    "${WORK}/src/generated/milca_vrmhd_dispatch.f90" RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "milca_vrmhd_dispatch.f90: regeneration differs")
endif()

execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${expected_directory}/milca_vrmhd_leaf_ids.f90"
    "${WORK}/src/generated/milca_vrmhd_leaf_ids.f90" RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "milca_vrmhd_leaf_ids.f90: regeneration differs")
endif()
