# Register kin6d_no_executable_stack over every KIN6D-owned executable and
# library. Dependencies are excluded: FetchContent subdirectories keep their
# own flags. Internal procedures passed as actual arguments still create
# trampolines; -ftrampoline-impl=heap (GCC >= 14) keeps them off the stack.
# Without it the stack is executable by construction, so the test exists
# only where the flag is available. The equilibrium CI tier builds through
# ci/cluster-toolchain.cmake with GCC 12.2: there the test is not registered
# and the binaries do have executable stacks. The guard covers GCC >= 14
# builds only.
function(_kin6d_stack_collect directory output)
    get_property(targets DIRECTORY "${directory}" PROPERTY BUILDSYSTEM_TARGETS)
    get_property(children DIRECTORY "${directory}" PROPERTY SUBDIRECTORIES)
    foreach(child IN LISTS children)
        cmake_path(IS_PREFIX PROJECT_SOURCE_DIR "${child}" NORMALIZE owned)
        cmake_path(IS_PREFIX PROJECT_BINARY_DIR "${child}" NORMALIZE generated)
        if(owned AND NOT generated)
            _kin6d_stack_collect("${child}" child_targets)
            list(APPEND targets ${child_targets})
        endif()
    endforeach()
    set(${output} ${targets} PARENT_SCOPE)
endfunction()

function(kin6d_stack_targets)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND NOT KIN6D_HAVE_HEAP_TRAMPOLINES)
        return()
    endif()
    _kin6d_stack_collect("${PROJECT_SOURCE_DIR}" targets)
    set(files)
    foreach(target IN LISTS targets)
        get_target_property(type ${target} TYPE)
        if(type MATCHES "^(EXECUTABLE|STATIC_LIBRARY|SHARED_LIBRARY|MODULE_LIBRARY)$")
            list(APPEND files "$<TARGET_FILE:${target}>")
        endif()
    endforeach()
    file(GENERATE OUTPUT "${CMAKE_BINARY_DIR}/stack-targets-$<CONFIG>.cmake"
        CONTENT "set(KIN6D_STACK_TARGETS \"${files}\")\n")
    find_program(KIN6D_READELF NAMES readelf llvm-readelf)
    add_test(NAME kin6d_no_executable_stack
        COMMAND ${CMAKE_COMMAND}
            "-DELF_FORMAT=${CMAKE_EXECUTABLE_FORMAT}"
            "-DREADELF=${KIN6D_READELF}"
            "-DTARGET_MANIFEST=${CMAKE_BINARY_DIR}/stack-targets-$<CONFIG>.cmake"
            -P "${PROJECT_SOURCE_DIR}/cmake/kin6d_no_executable_stack.cmake")
    set_tests_properties(kin6d_no_executable_stack PROPERTIES
        LABELS "core;fast" TIMEOUT 120
        SKIP_REGULAR_EXPRESSION "SKIP:.*require ELF")
endfunction()
