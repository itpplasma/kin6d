# Independent response readback/qualification; separate from coupled solver.
target_sources(kin6d PRIVATE src/response_calibration.f90
    src/generated/kin6d_response_calibration_leaves.f90)
add_executable(test_response_calibration test/test_response_calibration.f90)
target_link_libraries(test_response_calibration PRIVATE kin6d)
add_test(NAME response_spot_calibration_algebra COMMAND test_response_calibration)
add_executable(test_response_sheet_cylinder test/test_response_sheet_cylinder.f90)
target_link_libraries(test_response_sheet_cylinder PRIVATE kin6d)
add_test(NAME response_unit_sheet_spot_cylinder COMMAND test_response_sheet_cylinder)
set_tests_properties(response_unit_sheet_spot_cylinder PROPERTIES TIMEOUT 600
    ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
if(KIN6D_FORTSYM_BUILD)
    add_executable(gen_response_calibration gen/gen_response_calibration.f90)
    target_link_libraries(gen_response_calibration PRIVATE kin6d_gen_support)
endif()
target_sources(kin6d PRIVATE src/response_qualification.f90)
add_executable(test_response_qualification test/test_response_qualification.f90)
target_link_libraries(test_response_qualification PRIVATE kin6d)
add_test(NAME response_qualification_frozen_rules COMMAND test_response_qualification)
add_executable(response_qualify app/response_qualify.f90)
target_link_libraries(response_qualify PRIVATE kin6d kin6d_io)
add_executable(gs_response_calibration app/gs_response_calibration.f90)
target_link_libraries(gs_response_calibration PRIVATE kin6d kin6d_io)
if(Python3_EXECUTABLE)
    add_test(NAME response_qualification_external_adapter
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test/check_response_qualification_adapter.py
            $<TARGET_FILE:response_qualify>)
endif()
add_executable(test_response_strip_readback test/test_response_strip_readback.f90)
target_link_libraries(test_response_strip_readback PRIVATE kin6d)
add_test(NAME response_spot_resolved_single_ring COMMAND test_response_strip_readback)
set_tests_properties(response_spot_resolved_single_ring PROPERTIES TIMEOUT 300
    ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
configure_file(test/fixtures/response_calibration.nml.in
    ${CMAKE_CURRENT_BINARY_DIR}/response_calibration.nml @ONLY)
add_test(NAME response_calibration_deck_validation
    COMMAND gs_response_calibration --validate --config
        ${CMAKE_CURRENT_BINARY_DIR}/response_calibration.nml)
add_test(NAME response_calibration_deck_rejects_zero_mode
    COMMAND gs_response_calibration --validate --config
        ${CMAKE_CURRENT_SOURCE_DIR}/test/fixtures/response_calibration_invalid.nml)
set_tests_properties(response_calibration_deck_rejects_zero_mode PROPERTIES
    WILL_FAIL TRUE)
