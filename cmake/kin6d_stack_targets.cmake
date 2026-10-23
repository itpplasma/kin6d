# Called after all project targets exist. Dependencies retain their own flags:
# directory-wide options would also reach FetchContent subdirectories.
function(_kin6d_stack_collect directory output)
    get_property(targets DIRECTORY "${directory}" PROPERTY BUILDSYSTEM_TARGETS)
    get_property(children DIRECTORY "${directory}" PROPERTY SUBDIRECTORIES)
    foreach(child IN LISTS children)
        cmake_path(IS_PREFIX PROJECT_SOURCE_DIR "${child}" NORMALIZE owned)
        if(owned)
            _kin6d_stack_collect("${child}" child_targets)
            list(APPEND targets ${child_targets})
        endif()
    endforeach()
    set(${output} ${targets} PARENT_SCOPE)
endfunction()

function(kin6d_stack_targets)
    _kin6d_stack_collect("${PROJECT_SOURCE_DIR}" targets)
    set(files)
    set(build_targets)
    foreach(target IN LISTS targets)
        get_target_property(type ${target} TYPE)
        if(NOT type MATCHES "^(EXECUTABLE|STATIC_LIBRARY|SHARED_LIBRARY|MODULE_LIBRARY|OBJECT_LIBRARY)$")
            continue()
        endif()
        if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
            target_compile_options(${target} PRIVATE
                $<$<COMPILE_LANGUAGE:Fortran>:-Werror=trampolines>)
        endif()
        if(type STREQUAL "OBJECT_LIBRARY")
            list(APPEND files "$<TARGET_OBJECTS:${target}>")
        else()
            list(APPEND files "$<TARGET_FILE:${target}>")
        endif()
        list(APPEND build_targets ${target})
    endforeach()
    file(GENERATE OUTPUT "${CMAKE_BINARY_DIR}/stack-targets-$<CONFIG>.cmake"
        CONTENT "set(KIN6D_STACK_TARGETS \"${files}\")\n")
    find_program(KIN6D_READELF NAMES readelf llvm-readelf)
    add_custom_target(kin6d_check_no_executable_stack DEPENDS ${build_targets})
    add_test(NAME kin6d_no_executable_stack
        COMMAND ${CMAKE_COMMAND}
            "-DELF_FORMAT=${CMAKE_EXECUTABLE_FORMAT}"
            "-DREADELF=${KIN6D_READELF}"
            "-DTARGET_MANIFEST=${CMAKE_BINARY_DIR}/stack-targets-$<CONFIG>.cmake"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/kin6d_no_executable_stack.cmake")
    set_tests_properties(kin6d_no_executable_stack PROPERTIES
        LABELS "core;fast" SKIP_REGULAR_EXPRESSION "SKIP:.*require ELF")
endfunction()
