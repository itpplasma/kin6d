# Runs nn1_spec_case (nr=4, N=4, 64 planes) on the smoke grid, then checks it against direct Biot-Savart.
execute_process(COMMAND ${EXE} ${DATA}/case.nml ${DATA}/grid_smoke.txt ${OUT} 4 4 64
    RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
    message(FATAL_ERROR "nn1_spec_case failed: ${rc}")
endif()
execute_process(COMMAND ${PY} ${CHECK} ${OUT}.B.txt RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
    message(FATAL_ERROR "nn1_spec_case deviates from direct Biot-Savart")
endif()
