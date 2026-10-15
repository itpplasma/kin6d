file(MAKE_DIRECTORY "${WORK}")
file(WRITE "${WORK}/boundary.dat" "16
3.5 0
3.4619397663 0.1913417162
3.3535533906 0.3535533906
3.1913417162 0.4619397663
3 0.5
2.8086582838 0.4619397663
2.6464466094 0.3535533906
2.5380602337 0.1913417162
2.5 0
2.5380602337 -0.1913417162
2.6464466094 -0.3535533906
2.8086582838 -0.4619397663
3 -0.5
3.1913417162 -0.4619397663
3.3535533906 -0.3535533906
3.4619397663 -0.1913417162
")
file(WRITE "${WORK}/input.nml" "&gs_phase1
boundary_file='boundary.dat', boundary_nodes=16, max_area=0.08,
pprime=-1e5, ffprime=-1, f_edge=5, curved=.true., degree=2,
radial_points=8, angular_points=32, solve_wall_limit_s=20
/\n")
execute_process(COMMAND "${BINARY}" --config input.nml
    WORKING_DIRECTORY "${WORK}" RESULT_VARIABLE status
    OUTPUT_VARIABLE out ERROR_VARIABLE err TIMEOUT 25)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "production estimator run failed: ${out}${err}")
endif()
file(READ "${WORK}/summary.dat" summary)
if(NOT summary MATCHES "majorant_method equilibrated-rt" OR
   NOT summary MATCHES "majorant_claim practical estimate" OR
   NOT summary MATCHES "majorant_quadrature_status approximate" OR
   NOT summary MATCHES "majorant_scope solution-energy-error" OR
   NOT summary MATCHES "integrated_psi_goal_status approximate")
    message(FATAL_ERROR "missing truthful production error metadata: ${summary}")
endif()
find_program(NCDUMP ncdump REQUIRED)
execute_process(COMMAND "${NCDUMP}" -h "${WORK}/output/gs_phase1.nc"
    OUTPUT_VARIABLE header RESULT_VARIABLE status)
if(NOT status EQUAL 0 OR NOT header MATCHES "majorant_method = \"equilibrated-rt\"" OR
   NOT header MATCHES "majorant_quadrature_status = \"approximate\"" OR
   NOT header MATCHES "integrated_psi_radius" OR
   NOT header MATCHES "majorant:kin6d_status = \"practical-estimate\"")
    message(FATAL_ERROR "NetCDF lost estimator status or goal result: ${header}")
endif()
# The represented affine polygon and fixed source also admit an independent
# outward audit. This must change the reported claim only on the supported path.
file(READ "${WORK}/input.nml" input)
string(REPLACE "curved=.true." "curved=.false." input "${input}")
string(REPLACE "max_area=0.08" "max_area=0.02" input "${input}")
file(WRITE "${WORK}/input.nml" "${input}")
execute_process(COMMAND "${BINARY}" --config input.nml
    WORKING_DIRECTORY "${WORK}" RESULT_VARIABLE status
    OUTPUT_VARIABLE out ERROR_VARIABLE err TIMEOUT 25)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "production audit run failed: ${out}${err}")
endif()
file(READ "${WORK}/summary.dat" summary)
if(NOT summary MATCHES "majorant_method equilibrated-rt0-audit" OR
   NOT summary MATCHES "majorant_claim certified" OR
   NOT summary MATCHES "majorant_quadrature_status rigorous-rt0-audit" OR
   NOT summary MATCHES "majorant_certification_status certified-continuous")
    message(FATAL_ERROR "supported affine audit was not reported: ${summary}")
endif()
execute_process(COMMAND "${NCDUMP}" -h "${WORK}/output/gs_phase1.nc"
    OUTPUT_VARIABLE header RESULT_VARIABLE status)
if(NOT status EQUAL 0 OR
   NOT header MATCHES "majorant_certification_status = \"certified-continuous\"" OR
   NOT header MATCHES "majorant:kin6d_status = \"certified-continuous\"" OR
   NOT header MATCHES "psi:kin6d_status = \"candidate\"")
    message(FATAL_ERROR "NetCDF omitted continuous audit status: ${header}")
endif()
