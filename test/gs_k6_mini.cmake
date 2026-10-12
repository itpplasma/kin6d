# K6-mini regression (docs/equilibrium-core.md, migration M0): gs_phase1 on
# the K6 law and path with a coarse mesh, then a native field query, compared
# with the fixture to 1e-12 relative by CHECK. The fixture was produced by
# this script at the M0 base (origin/main 260ce1d) with OMP_NUM_THREADS=4;
# boundary.dat is every eighth sample of the K6 boundary, source.dat the K6
# cubic-source law unchanged.
file(REMOVE_RECURSE "${WORK}")
file(MAKE_DIRECTORY "${WORK}")
foreach(name boundary.dat source.dat input.nml points.dat)
    file(COPY "${FIXTURES}/${name}" DESTINATION "${WORK}")
endforeach()
execute_process(COMMAND "${APP}" --config input.nml WORKING_DIRECTORY "${WORK}"
    OUTPUT_FILE "${WORK}/native.log" ERROR_FILE "${WORK}/native.err"
    RESULT_VARIABLE status)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "K6-mini gs_phase1 failed (${status}); see ${WORK}/native.log")
endif()
execute_process(COMMAND "${APP}" --evaluate mesh.dat 1 0 profile-law.dat
    WORKING_DIRECTORY "${WORK}" INPUT_FILE "${WORK}/points.dat"
    OUTPUT_FILE "${WORK}/field.dat" RESULT_VARIABLE status)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "K6-mini field query failed (${status})")
endif()
execute_process(COMMAND "${CHECK}" compare "${WORK}" "${FIXTURES}"
    RESULT_VARIABLE status)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "K6-mini outputs differ from the fixture")
endif()
