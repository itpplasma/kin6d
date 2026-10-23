# Optional research spike: generic tensor machinery is compiled from FortFEM,
# never copied into KIN6D. The default triangle solver is untouched.
add_library(kin6d_gll_mf STATIC
    ${fortfem_SOURCE_DIR}/src/discretisation/tensor_line_tables.f90
    ${fortfem_SOURCE_DIR}/src/discretisation/gll_line.f90
    ${fortfem_SOURCE_DIR}/src/discretisation/tensor_mesh.f90
    ${fortfem_SOURCE_DIR}/src/discretisation/tensor_contractions.f90
    ${fortfem_SOURCE_DIR}/src/discretisation/line_mesh.f90
    ${fortfem_SOURCE_DIR}/src/generated/fortfem_tensor_metric_products.f90
    ${fortfem_SOURCE_DIR}/src/discretisation/tensor_operator.f90
    ${fortfem_SOURCE_DIR}/src/discretisation/tensor_solver.f90)
target_link_libraries(kin6d_gll_mf PUBLIC kin6d_sparse_direct OpenMP::OpenMP_Fortran)
add_executable(benchmark_ideal_gll
    ${PROJECT_SOURCE_DIR}/src/generated/kin6d_ideal_gll_leaves.f90
    ideal_gll.f90 benchmark_gll.f90)
target_link_libraries(benchmark_ideal_gll PRIVATE kin6d_mhd kin6d_gll_mf)
set_source_files_properties(${PROJECT_SOURCE_DIR}/src/generated/kin6d_ideal_gll_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
