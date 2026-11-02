# Layered KIN6D libraries (ROADMAP "Build structure and fast development loop").
#
#   L0  kin6d_core         intervals, ODE, species, sort keys, statistics, laboratory
#   L1  kin6d_fem2d        FortFEM slice, triangle bases, CSR/direct solvers, GS options
#   L1  kin6d_pairs        Ewald, PME, N-body, pair laws, radial DVR, collision kernels
#   L2  kin6d_equilibrium  GS, equilibrium, 3D coils, equilibrated bounds
#   L2  kin6d_mhd          MHD, ideal maps, response, vacuum surfaces, nonnested
#   L2  kin6d_classical, kin6d_quantum, kin6d_continuum
#   L3  applications; kin6d_io (FortIO results, animation export) and plotting
#
# A layer uses only layers below it; among L2 classes only mhd uses equilibrium.
# Sources are assigned with target_sources(<layer> ...) in CMakeLists.txt; the
# kin6d_layering test checks the Fortran use graph against that assignment.
# `kin6d` is the umbrella INTERFACE target for consumers that want everything.

set(KIN6D_LAYER_CLASSES core fem2d pairs equilibrium mhd classical quantum continuum)
# Test labels: the layer classes plus apps for application-only tests.
set(KIN6D_TEST_CLASSES ${KIN6D_LAYER_CLASSES} apps)
set(KIN6D_LAYER_LIBRARIES)
foreach(layer IN LISTS KIN6D_LAYER_CLASSES)
    add_library(kin6d_${layer})
    list(APPEND KIN6D_LAYER_LIBRARIES kin6d_${layer})
endforeach()
target_link_libraries(kin6d_core PUBLIC kin6d_intervals kin6d_ode)
target_link_libraries(kin6d_fem2d PUBLIC kin6d_core)
target_link_libraries(kin6d_pairs PUBLIC kin6d_core)
target_link_libraries(kin6d_equilibrium PUBLIC kin6d_fem2d)
target_link_libraries(kin6d_mhd PUBLIC kin6d_equilibrium kin6d_fem2d)
target_link_libraries(kin6d_classical PUBLIC kin6d_pairs)
target_link_libraries(kin6d_quantum PUBLIC kin6d_pairs)
target_link_libraries(kin6d_continuum PUBLIC kin6d_core)
add_library(kin6d INTERFACE)
target_link_libraries(kin6d INTERFACE ${KIN6D_LAYER_LIBRARIES})

option(KIN6D_BUILD_EXTRAS
    "Build applications, plotting and result IO with the default target" ON)

# Collect "module name -> layer" for the sources of the given targets.
function(_kin6d_layer_scan_modules layer)
    foreach(target IN LISTS ARGN)
        if(NOT TARGET ${target})
            continue()
        endif()
        get_target_property(sources ${target} SOURCES)
        get_target_property(source_dir ${target} SOURCE_DIR)
        foreach(source IN LISTS sources)
            cmake_path(ABSOLUTE_PATH source BASE_DIRECTORY "${source_dir}")
            if(NOT source MATCHES "\\.[fF]90$" OR NOT EXISTS "${source}")
                continue()
            endif()
            set_property(GLOBAL APPEND PROPERTY KIN6D_LAYER_SOURCES "${layer}=${source}")
            file(STRINGS "${source}" lines REGEX
                "^[ \t]*[Mm][Oo][Dd][Uu][Ll][Ee][ \t]+[A-Za-z0-9_]+[ \t]*(!.*)?$")
            foreach(line IN LISTS lines)
                string(REGEX REPLACE "^[ \t]*[A-Za-z]+[ \t]+([A-Za-z0-9_]+).*$" "\\1"
                    module "${line}")
                string(TOLOWER "${module}" module)
                if(NOT module STREQUAL "procedure")
                    set_property(GLOBAL PROPERTY KIN6D_MODULE_LAYER_${module} ${layer})
                endif()
            endforeach()
        endforeach()
    endforeach()
endfunction()

# Layer targets whose modules a target's own sources use.
function(_kin6d_layer_needs target out)
    get_target_property(sources ${target} SOURCES)
    get_target_property(source_dir ${target} SOURCE_DIR)
    set(needs)
    foreach(source IN LISTS sources)
        cmake_path(ABSOLUTE_PATH source BASE_DIRECTORY "${source_dir}")
        if(NOT source MATCHES "\\.[fF]90$" OR NOT EXISTS "${source}")
            continue()
        endif()
        file(STRINGS "${source}" lines REGEX
            "^[ \t]*([Uu][Ss][Ee]|[Ss][Uu][Bb][Mm][Oo][Dd][Uu][Ll][Ee])[ \t,:(]")
        foreach(line IN LISTS lines)
            string(REGEX REPLACE
                "^[ \t]*[A-Za-z]+[ \t]*(,[ \t]*[A-Za-z_]+[ \t]*)?(::)?[ \t(]*([A-Za-z0-9_]+).*$"
                "\\3" module "${line}")
            string(TOLOWER "${module}" module)
            get_property(layer GLOBAL PROPERTY KIN6D_MODULE_LAYER_${module})
            if(layer AND NOT layer STREQUAL "apps")
                list(APPEND needs kin6d_${layer})
            endif()
        endforeach()
    endforeach()
    list(REMOVE_DUPLICATES needs)
    set(${out} "${needs}" PARENT_SCOPE)
endfunction()

# Replace a link to the umbrella by the layers the target's sources use, so a
# test or application rebuilds and relinks only against what it needs.
# A target that reached the umbrella only through a library it links (a test
# fixture) gets the layers its own sources use added to its links.
function(_kin6d_layer_narrow target carriers)
    _kin6d_layer_needs(${target} needs)
    get_target_property(libraries ${target} LINK_LIBRARIES)
    if(libraries AND NOT "kin6d" IN_LIST libraries)
        foreach(library IN LISTS libraries)
            if(library IN_LIST carriers)
                set_property(TARGET ${target} APPEND PROPERTY LINK_LIBRARIES ${needs})
                break()
            endif()
        endforeach()
    endif()
    foreach(property LINK_LIBRARIES INTERFACE_LINK_LIBRARIES)
        get_target_property(libraries ${target} ${property})
        if(NOT libraries)
            continue()
        endif()
        set(narrowed)
        set(changed FALSE)
        foreach(library IN LISTS libraries)
            if(library STREQUAL "kin6d")
                list(APPEND narrowed ${needs})
                set(changed TRUE)
            elseif(library STREQUAL "$<LINK_ONLY:kin6d>")
                foreach(need IN LISTS needs)
                    list(APPEND narrowed "$<LINK_ONLY:${need}>")
                endforeach()
                set(changed TRUE)
            else()
                list(APPEND narrowed "${library}")
            endif()
        endforeach()
        if(changed)
            set_property(TARGET ${target} PROPERTY ${property} "${narrowed}")
        endif()
    endforeach()
endfunction()

# Label one test with its model class(es) and tier and record the executables
# its kin6d_check_<class> target builds. Tests not configured here are skipped.
function(kin6d_test_labels test tier)
    cmake_parse_arguments(PARSE_ARGV 2 arg "" "" "CLASSES;TARGETS")
    if(NOT TEST ${test})
        return()
    endif()
    set_property(TEST ${test} APPEND PROPERTY LABELS ${arg_CLASSES} ${tier})
    foreach(target IN LISTS arg_TARGETS)
        if(NOT TARGET ${target})
            continue()
        endif()
        get_target_property(excluded ${target} EXCLUDE_FROM_ALL)
        if(excluded)
            set_property(TEST ${test} PROPERTY DISABLED TRUE)
            return()
        endif()
    endforeach()
    foreach(class IN LISTS arg_CLASSES)
        set_property(GLOBAL APPEND PROPERTY KIN6D_CHECK_TARGETS_${class} ${arg_TARGETS})
    endforeach()
endfunction()

# Called once at the end of CMakeLists.txt, after every target and test exists.
function(kin6d_layers_finalize)
    foreach(layer IN LISTS KIN6D_LAYER_CLASSES)
        _kin6d_layer_scan_modules(${layer} kin6d_${layer})
    endforeach()
    _kin6d_layer_scan_modules(core kin6d_intervals kin6d_ode kin6d_fft
        kin6d_libneo_coils)
    _kin6d_layer_scan_modules(fem2d kin6d_spatial kin6d_sparse_direct kin6d_feec_n
        fortsparse)
    _kin6d_layer_scan_modules(apps kin6d_io kin6d_geqdsk fortplot)
    get_directory_property(targets DIRECTORY ${CMAKE_SOURCE_DIR} BUILDSYSTEM_TARGETS)
    set(narrow)
    set(carriers kin6d)
    foreach(target IN LISTS targets)
        get_target_property(type ${target} TYPE)
        if((type STREQUAL "EXECUTABLE" OR type STREQUAL "STATIC_LIBRARY")
                AND NOT target IN_LIST KIN6D_LAYER_LIBRARIES)
            list(APPEND narrow ${target})
        endif()
    endforeach()
    # Libraries that pass the umbrella on to their consumers, transitively.
    set(grown TRUE)
    while(grown)
        set(grown FALSE)
        foreach(target IN LISTS narrow)
            get_target_property(interface ${target} INTERFACE_LINK_LIBRARIES)
            if(NOT target IN_LIST carriers AND interface)
                foreach(library IN LISTS interface)
                    if(library IN_LIST carriers)
                        list(APPEND carriers ${target})
                        set(grown TRUE)
                        break()
                    endif()
                endforeach()
            endif()
        endforeach()
    endwhile()
    foreach(target IN LISTS narrow)
        _kin6d_layer_narrow(${target} "${carriers}")
    endforeach()
    if(NOT KIN6D_BUILD_EXTRAS)
        # Applications, plotting and result IO leave the default build; their
        # tests are disabled by kin6d_test_labels.
        foreach(target IN LISTS targets)
            get_target_property(type ${target} TYPE)
            get_target_property(libraries ${target} LINK_LIBRARIES)
            get_target_property(sources ${target} SOURCES)
            if(NOT type MATCHES "EXECUTABLE|STATIC_LIBRARY")
                continue()
            endif()
            if(target STREQUAL "kin6d_io" OR "${sources}" MATCHES "(^|;)app/"
                    OR "${libraries}" MATCHES "(^|;)(fortplot|kin6d_io|kin6d_demo_style)(;|$)")
                set_property(TARGET ${target} PROPERTY EXCLUDE_FROM_ALL TRUE)
            endif()
        endforeach()
        foreach(dependency fortplot fortio)
            if(DEFINED ${dependency}_BINARY_DIR AND IS_DIRECTORY "${${dependency}_BINARY_DIR}")
                set_property(DIRECTORY "${${dependency}_SOURCE_DIR}"
                    PROPERTY EXCLUDE_FROM_ALL TRUE)
            endif()
        endforeach()
    endif()
    include(${CMAKE_SOURCE_DIR}/cmake/kin6d_test_labels.cmake)
    foreach(class IN LISTS KIN6D_TEST_CLASSES)
        get_property(check_targets GLOBAL PROPERTY KIN6D_CHECK_TARGETS_${class})
        set(existing)
        foreach(target IN LISTS check_targets)
            if(TARGET ${target})
                list(APPEND existing ${target})
            endif()
        endforeach()
        list(REMOVE_DUPLICATES existing)
        add_custom_target(kin6d_check_${class} DEPENDS ${existing}
            COMMENT "Built the ${class} test executables; run: ctest -L ${class} -LE slow")
    endforeach()
    # Partition for the layering test: one "layer=source" line per source.
    get_property(layer_sources GLOBAL PROPERTY KIN6D_LAYER_SOURCES)
    string(REPLACE ";" "\n" layer_sources "${layer_sources}")
    file(WRITE ${CMAKE_BINARY_DIR}/kin6d_layer_sources.txt "${layer_sources}\n")
    add_test(NAME kin6d_layering
        COMMAND ${CMAKE_COMMAND}
            -DPARTITION=${CMAKE_BINARY_DIR}/kin6d_layer_sources.txt
            -DCTEST=${CMAKE_CTEST_COMMAND} -DBUILD_DIR=${CMAKE_BINARY_DIR}
            "-DCLASSES=${KIN6D_TEST_CLASSES}"
            -P ${CMAKE_SOURCE_DIR}/cmake/kin6d_check_layers.cmake)
    set_tests_properties(kin6d_layering PROPERTIES LABELS "core;fast")
endfunction()
