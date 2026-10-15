file(REMOVE_RECURSE "${WORK}")
file(MAKE_DIRECTORY "${WORK}")
file(WRITE "${WORK}/input.nml"
    "&grad_shafranov solve_mode='prescribed_constant_q',element_order=2,r0=4,a=.4,b0=.001,q_seed=0,max_area=.0025,boundary_nodes=63,inverse_flux_points=5,inverse_angular_points=128,q_target=1.5,toroidal_flux_target=1,q_tolerance=.005,nonlinear_tolerance=1e-6 /\n")
execute_process(COMMAND "${APP}" --config input.nml WORKING_DIRECTORY "${WORK}"
    RESULT_VARIABLE status OUTPUT_VARIABLE out ERROR_VARIABLE err TIMEOUT 600)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "native inverse app failed: ${out}${err}")
endif()
if(NOT out MATCHES "prescribed q/Phi discrete indicators")
    message(FATAL_ERROR "app did not execute prescribed-q entry point")
endif()
find_program(NCDUMP ncdump REQUIRED)
execute_process(COMMAND "${NCDUMP}" -v axis_q,q_target,Phi_target,F_edge,residual
    "${WORK}/output/grad_shafranov.nc" RESULT_VARIABLE status OUTPUT_VARIABLE data)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "independent NetCDF reader failed")
endif()
foreach(name axis_q q_target Phi_target F_edge residual)
    string(REGEX MATCH "[\n ]${name} = ([^ ;]+)" match "${data}")
    if(NOT match)
        message(FATAL_ERROR "missing native inverse result ${name}")
    endif()
    set(${name} "${CMAKE_MATCH_1}")
endforeach()
if(NOT axis_q GREATER 1.49 OR NOT axis_q LESS 1.51 OR NOT q_target EQUAL 1.5)
    message(FATAL_ERROR "native app prescribed safety factor differs from independent axis oracle")
endif()
if(NOT Phi_target EQUAL 1 OR NOT residual LESS 1e-6)
    message(FATAL_ERROR "native app dropped Phi or weak GS stopping")
endif()
if(NOT F_edge GREATER 1)
    message(FATAL_ERROR "inverse used the forward b0/source seed instead of solving F")
endif()
# The inverse state is estimated through the fixed-source equilibrated path;
# sampled quadrature keeps the claim practical, never certified.
execute_process(COMMAND "${NCDUMP}" -h "${WORK}/output/grad_shafranov.nc"
    RESULT_VARIABLE status OUTPUT_VARIABLE header)
if(NOT status EQUAL 0 OR NOT header MATCHES ":majorant_claim = \"practical estimate\""
        OR NOT header MATCHES ":majorant_method = \"equilibrated-rt\""
        OR NOT header MATCHES ":majorant_quadrature_status = \"approximate\"")
    message(FATAL_ERROR "inverse result lacks the practical fixed-source estimate claim")
endif()
execute_process(COMMAND "${NCDUMP}" -v majorant,majorant_stability,majorant_bound,majorant_status
    "${WORK}/output/grad_shafranov.nc" RESULT_VARIABLE status OUTPUT_VARIABLE data)
foreach(name majorant majorant_stability majorant_bound majorant_status)
    string(REGEX MATCH "[\n ]${name} = ([^ ;]+)" match "${data}")
    if(NOT match)
        message(FATAL_ERROR "missing inverse estimate ${name}")
    endif()
    set(${name} "${CMAKE_MATCH_1}")
endforeach()
if(NOT majorant_status EQUAL 5 OR NOT majorant GREATER 0 OR NOT majorant_stability GREATER 0
        OR NOT majorant_bound GREATER 0)
    message(FATAL_ERROR "inverse fixed-source estimate is not a positive practical estimate")
endif()
file(WRITE "${WORK}/unsupported.nml"
    "&grad_shafranov solve_mode='prescribed_constant_q',element_order=2,pressure_mode='finite_beta',toroidal_flux_target=1 /\n")
execute_process(COMMAND "${APP}" --config unsupported.nml WORKING_DIRECTORY "${WORK}"
    RESULT_VARIABLE status OUTPUT_VARIABLE out ERROR_VARIABLE err TIMEOUT 120)
if(status EQUAL 0 OR NOT "${out}${err}" MATCHES "requires zero pressure")
    message(FATAL_ERROR "unsupported finite-pressure inverse silently accepted")
endif()
