# Cap the output CTest records per test (decision 54). CTest writes the
# complete output of every test into Testing/Temporary/LastTest.log, and the
# CTEST_CUSTOM_MAXIMUM_*_TEST_OUTPUT_SIZE limits only shorten dashboard
# submissions, not that file. Every executable target therefore runs its tests
# through cmake/kin6d_test_output_cap.sh, which keeps the first
# KIN6D_TEST_OUTPUT_LIMIT bytes and the last 200 lines and preserves the exit
# status. Include this file after the last add_executable of the project.
set(KIN6D_TEST_OUTPUT_LIMIT 16777216 CACHE STRING
    "Bytes of output kept per test before truncation (0 disables the cap)")
if(BUILD_TESTING AND KIN6D_TEST_OUTPUT_LIMIT GREATER 0)
    find_program(KIN6D_BASH bash REQUIRED)
    set(_kin6d_cap_launcher "${KIN6D_BASH}"
        "${CMAKE_CURRENT_LIST_DIR}/kin6d_test_output_cap.sh"
        "${KIN6D_TEST_OUTPUT_LIMIT}")
    get_property(_kin6d_cap_targets DIRECTORY PROPERTY BUILDSYSTEM_TARGETS)
    foreach(_t IN LISTS _kin6d_cap_targets)
        get_target_property(_type ${_t} TYPE)
        if(NOT _type STREQUAL "EXECUTABLE")
            continue()
        endif()
        if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.29)
            set_property(TARGET ${_t} PROPERTY TEST_LAUNCHER ${_kin6d_cap_launcher})
        else()
            # Before TEST_LAUNCHER, add_test prefixes the emulator of a target.
            # No add_custom_command of this project runs a built executable.
            set_property(TARGET ${_t} PROPERTY CROSSCOMPILING_EMULATOR
                ${_kin6d_cap_launcher})
        endif()
    endforeach()
endif()
