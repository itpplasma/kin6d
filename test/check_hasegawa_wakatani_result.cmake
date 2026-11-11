# Output contract of the drift-turbulence branch, run by CTest with -P.
# Inputs: APP (hasegawa_wakatani), CHECK (test_turbulence_output), WORK.
# The application runs a small saturated case; the checker reads the result
# file through FortIO's HDF5 reader.
file(REMOVE_RECURSE "${WORK}")
file(MAKE_DIRECTORY "${WORK}")
file(WRITE "${WORK}/case.nml" "&hasegawa_wakatani
 alpha = 0.5, kappa = 1.0, d = 2.0e-3, order = 2, modified = .false.,
 nx = 32, ny = 32, k0 = 0.3, dt = 0.025, time_end = 500.0, average_from = 100.0,
 sample_steps = 16, segment = 256, modes = 10, amplitude = 1.0e-2, seed = 3,
 rho_s_cm = 0.2, inverse_omega_ci_s = 1.0e-8, te_over_e_statvolt = 2.5 /
")
foreach(step APP CHECK)
    if(step STREQUAL "APP")
        set(command "${APP}" --config case.nml)
    else()
        set(command "${CHECK}")
    endif()
    execute_process(COMMAND ${command} WORKING_DIRECTORY "${WORK}"
        RESULT_VARIABLE code OUTPUT_VARIABLE out ERROR_VARIABLE err)
    message(STATUS "${out}")
    if(NOT code EQUAL 0)
        message(FATAL_ERROR "${command} failed:\n${out}${err}")
    endif()
endforeach()
