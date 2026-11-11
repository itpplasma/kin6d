# Drift-wave turbulence branch (docs/drift-turbulence.md).
target_sources(milca_continuum PRIVATE src/generated/milca_drift_fluid_leaves.f90
    src/continuum/drift_dispersion.f90 src/continuum/drift_spectral.f90
    src/continuum/hasegawa_wakatani.f90 src/continuum/turbulence_statistics.f90)
set_source_files_properties(src/generated/milca_drift_fluid_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
if(MILCA_FORTSYM_BUILD)
    add_executable(gen_drift_fluid gen/gen_drift_fluid.f90)
    target_link_libraries(gen_drift_fluid PRIVATE milca_gen_support)
    add_test(NAME drift_fluid_leaf_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_drift_fluid>
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/milca_drift_fluid_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/drift-fluid-regeneration
        -DLEAF=milca_drift_fluid_leaves.f90
        -P ${CMAKE_SOURCE_DIR}/cmake/milca_check_leaf_regeneration.cmake)
    milca_test_labels(drift_fluid_leaf_regeneration fast CLASSES continuum
        TARGETS gen_drift_fluid)
endif()

add_executable(test_drift_dispersion test/test_drift_dispersion.f90)
target_link_libraries(test_drift_dispersion PRIVATE milca_continuum)
add_test(NAME drift_dispersion_fluid_roots COMMAND test_drift_dispersion)
set_tests_properties(drift_dispersion_fluid_roots PROPERTIES TIMEOUT 60)
milca_test_labels(drift_dispersion_fluid_roots fast CLASSES continuum
    TARGETS test_drift_dispersion)

add_executable(test_drift_spectral test/test_drift_spectral.f90)
target_link_libraries(test_drift_spectral PRIVATE milca_continuum)
add_test(NAME drift_spectral_bracket COMMAND test_drift_spectral)
set_tests_properties(drift_spectral_bracket PROPERTIES TIMEOUT 60)
milca_test_labels(drift_spectral_bracket fast CLASSES continuum
    TARGETS test_drift_spectral)

add_executable(test_hasegawa_wakatani test/test_hasegawa_wakatani.f90)
target_link_libraries(test_hasegawa_wakatani PRIVATE milca_continuum)
add_test(NAME hasegawa_wakatani_oracles COMMAND test_hasegawa_wakatani)
set_tests_properties(hasegawa_wakatani_oracles PROPERTIES TIMEOUT 120)
milca_test_labels(hasegawa_wakatani_oracles fast CLASSES continuum
    TARGETS test_hasegawa_wakatani)

add_executable(test_turbulence_statistics test/test_turbulence_statistics.f90)
target_link_libraries(test_turbulence_statistics PRIVATE milca_continuum)
add_test(NAME turbulence_statistics_oracles COMMAND test_turbulence_statistics)
set_tests_properties(turbulence_statistics_oracles PROPERTIES TIMEOUT 120)
milca_test_labels(turbulence_statistics_oracles fast CLASSES continuum
    TARGETS test_turbulence_statistics)

target_sources(milca_io PRIVATE src/apps/turbulence_output.f90)
add_executable(hasegawa_wakatani app/hasegawa_wakatani.f90)
target_link_libraries(hasegawa_wakatani PRIVATE milca_io)
add_executable(test_turbulence_output test/test_turbulence_output.f90)
target_link_libraries(test_turbulence_output PRIVATE milca_io)
add_test(NAME hasegawa_wakatani_result_contract COMMAND ${CMAKE_COMMAND}
    -DAPP=$<TARGET_FILE:hasegawa_wakatani> -DCHECK=$<TARGET_FILE:test_turbulence_output>
    -DWORK=${CMAKE_BINARY_DIR}/hasegawa-wakatani-result
    -P ${CMAKE_SOURCE_DIR}/test/check_hasegawa_wakatani_result.cmake)
set_tests_properties(hasegawa_wakatani_result_contract PROPERTIES TIMEOUT 300
    ENVIRONMENT "OMP_NUM_THREADS=1")
milca_test_labels(hasegawa_wakatani_result_contract fast CLASSES continuum apps
    TARGETS hasegawa_wakatani test_turbulence_output)

# Saturated flux against published values (benchmarks/drift-turbulence). The
# modified model, whose zonal flows suppress the flux, has to be rejected.
add_executable(test_hasegawa_wakatani_benchmark test/test_hasegawa_wakatani_benchmark.f90)
target_link_libraries(test_hasegawa_wakatani_benchmark PRIVATE milca_io)
foreach(variant flux zonal_mutant)
    if(variant STREQUAL "flux")
        set(modified .false.)
        set(reject FALSE)
    else()
        set(modified .true.)
        set(reject TRUE)
    endif()
    add_test(NAME hasegawa_wakatani_benchmark_${variant} COMMAND ${CMAKE_COMMAND}
        -DAPP=$<TARGET_FILE:hasegawa_wakatani>
        -DCHECK=$<TARGET_FILE:test_hasegawa_wakatani_benchmark>
        -DWORK=${CMAKE_BINARY_DIR}/hasegawa-wakatani-benchmark-${variant}
        -DMODIFIED=${modified} -DREJECT=${reject}
        -P ${CMAKE_SOURCE_DIR}/test/check_hasegawa_wakatani_benchmark.cmake)
    set_tests_properties(hasegawa_wakatani_benchmark_${variant} PROPERTIES TIMEOUT 3600
        ENVIRONMENT "OMP_NUM_THREADS=1")
    milca_test_labels(hasegawa_wakatani_benchmark_${variant} slow CLASSES continuum apps
        TARGETS hasegawa_wakatani test_hasegawa_wakatani_benchmark)
endforeach()
