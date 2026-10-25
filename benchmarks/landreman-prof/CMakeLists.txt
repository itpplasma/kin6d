# Measurement-only target (perf/landreman-w3-20261025), never part of a PR.
add_executable(landreman_prof EXCLUDE_FROM_ALL
    ${CMAKE_CURRENT_LIST_DIR}/kin6d_map_3d_readback_leaves.f90
    ${CMAKE_CURRENT_LIST_DIR}/map_3d_equilibrium.f90
    ${CMAKE_CURRENT_LIST_DIR}/equilibrium_3d_driver.f90
    ${CMAKE_CURRENT_LIST_DIR}/map_spectrum_probe.f90
    ${CMAKE_CURRENT_LIST_DIR}/landreman_prof_main.F90)
target_link_libraries(landreman_prof PRIVATE kin6d_equilibrium kin6d)
if(KIN6D_PROF_KRYLOV)
    target_compile_definitions(landreman_prof PRIVATE KIN6D_PROF_KRYLOV)
endif()
