# Saturated Hasegawa-Wakatani benchmark (benchmarks/drift-turbulence), run by
# CTest with -P. Inputs: APP (hasegawa_wakatani), CHECK
# (test_hasegawa_wakatani_benchmark), WORK, ALPHA, MODIFIED (.false. or .true.)
# and REJECT (TRUE when the acceptance check has to fail).
file(REMOVE_RECURSE "${WORK}")
file(MAKE_DIRECTORY "${WORK}")
file(WRITE "${WORK}/case.nml" "&hasegawa_wakatani
 alpha = ${ALPHA}, kappa = 1.0, d = 1.0e-5, order = 3, modified = ${MODIFIED},
 nx = 128, ny = 128, k0 = 0.15, dt = 0.02, time_end = 1000.0, average_from = 300.0,
 sample_steps = 10, segment = 512, modes = 16, amplitude = 1.0e-2, seed = 1 /
")
execute_process(COMMAND "${APP}" --config case.nml WORKING_DIRECTORY "${WORK}"
    RESULT_VARIABLE code OUTPUT_VARIABLE out ERROR_VARIABLE err)
message(STATUS "${out}")
if(NOT code EQUAL 0)
    message(FATAL_ERROR "${APP} failed:\n${out}${err}")
endif()
execute_process(COMMAND "${CHECK}" WORKING_DIRECTORY "${WORK}"
    RESULT_VARIABLE code OUTPUT_VARIABLE out ERROR_VARIABLE err)
message(STATUS "${out}")
if(REJECT AND code EQUAL 0)
    message(FATAL_ERROR "the acceptance check passed a run that it has to reject")
elseif(NOT REJECT AND NOT code EQUAL 0)
    message(FATAL_ERROR "${CHECK} failed:\n${out}${err}")
endif()
