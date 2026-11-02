# Ideal current sheets at rational surfaces (docs/current-sheet-equilibrium.md).
target_sources(kin6d_mhd PRIVATE src/generated/kin6d_current_sheet_leaves.f90
    src/mhd/current_sheet_oracles.f90)
if(KIN6D_FORTSYM_BUILD)
    add_executable(gen_current_sheet gen/gen_current_sheet.f90)
    target_link_libraries(gen_current_sheet PRIVATE kin6d_gen_support)
    add_test(NAME current_sheet_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_current_sheet>
        -DFAMILY=current_sheet
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/kin6d_current_sheet_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/current-sheet-regeneration
        -P ${CMAKE_SOURCE_DIR}/cmake/kin6d_check_e21_regeneration.cmake)
    set_tests_properties(current_sheet_regeneration PROPERTIES TIMEOUT 60)
    kin6d_test_labels(current_sheet_regeneration fast CLASSES mhd
        TARGETS gen_current_sheet)
endif()

add_executable(test_current_sheet_oracles test/test_current_sheet_oracles.f90)
target_link_libraries(test_current_sheet_oracles PRIVATE kin6d_mhd)
foreach(case slab_identities slab_fe screw_pinch)
    add_test(NAME current_sheet_oracle_${case}
        COMMAND test_current_sheet_oracles ${case})
    set_tests_properties(current_sheet_oracle_${case} PROPERTIES TIMEOUT 60)
    kin6d_test_labels(current_sheet_oracle_${case} fast CLASSES mhd
        TARGETS test_current_sheet_oracles)
endforeach()
