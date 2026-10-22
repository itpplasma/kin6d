cmake_minimum_required(VERSION 3.20)
# IH2 pilot (optional accelerator): nested-solution charts on FortFEM's MF1
# matrix-free tensor operator. Standalone so KIN6D's main FortFEM pin is not
# changed; the lab pin moves to FortFEM main once lazy-fortran/fortfem#110 lands.
project(kin6d_c1_charts LANGUAGES Fortran C CXX)
include(FetchContent)
set(KIN6D_C1_FORTFEM_REVISION "52cd3be29c8b3215441ed716a9a0f8cafcfe822c"
    CACHE STRING "FortFEM lab/mf1 revision (PR #110 head)")
set(FORTFEM_NATIVE_SPARSE_ONLY ON CACHE BOOL "" FORCE)
FetchContent_Declare(fortfem
    GIT_REPOSITORY https://github.com/lazy-fortran/fortfem.git
    GIT_TAG ${KIN6D_C1_FORTFEM_REVISION})
FetchContent_MakeAvailable(fortfem)
find_package(OpenMP REQUIRED COMPONENTS Fortran)

add_library(c1_nested_chart STATIC nested_chart.f90)
set_target_properties(c1_nested_chart PROPERTIES
    Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/chart_modules")
target_include_directories(c1_nested_chart PUBLIC "${PROJECT_BINARY_DIR}/chart_modules")

add_executable(c1_charts c1_charts.f90)
target_link_libraries(c1_charts PRIVATE c1_nested_chart fortfem_tensor_solver
    OpenMP::OpenMP_Fortran)

enable_testing()
add_executable(test_nested_chart test_nested_chart.f90)
target_link_libraries(test_nested_chart PRIVATE c1_nested_chart)
add_test(NAME c1_nested_chart COMMAND test_nested_chart)
set_tests_properties(c1_nested_chart PROPERTIES LABELS "c1;fast" TIMEOUT 60)
set(KIN6D_C1_MAP "" CACHE FILEPATH "Optional fitted nested map for the file test")
if(KIN6D_C1_MAP)
    add_test(NAME c1_nested_chart_file COMMAND test_nested_chart ${KIN6D_C1_MAP})
    add_test(NAME c1_charts_smoke
        COMMAND c1_charts ${KIN6D_C1_MAP} aligned 2 2 2 1e4 1e-10)
    set_tests_properties(c1_nested_chart_file c1_charts_smoke PROPERTIES
        LABELS "c1" TIMEOUT 300 ENVIRONMENT "OMP_NUM_THREADS=4")
endif()
