file(MAKE_DIRECTORY "${WORK}")
file(WRITE "${WORK}/mesh.dat" "# kin6d-gs-mesh-v2: P2 SI R Z psi boundary
6 1
1 0 1 1
2 0 4 1
1 1 2 1
1.5 0 2.25 1
1.5 0.5 2.5 1
1 0.5 1.25 1
# connectivity
1 2 3 4 5 6
")
set(request "3\n1.2 0.1\n1 0\n3 0\n")
file(WRITE "${WORK}/request.txt" "${request}")
execute_process(COMMAND "${BINARY}" --evaluate mesh.dat 5 0
    WORKING_DIRECTORY "${WORK}" INPUT_FILE "${WORK}/request.txt"
    OUTPUT_VARIABLE single RESULT_VARIABLE status)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "single field query failed")
endif()
file(WRITE "${WORK}/request.txt" "${request}0\n${request}")
execute_process(COMMAND "${BINARY}" --evaluate-stream mesh.dat 5 0
    WORKING_DIRECTORY "${WORK}" INPUT_FILE "${WORK}/request.txt"
    OUTPUT_VARIABLE repeated RESULT_VARIABLE status TIMEOUT 10)
if(NOT status EQUAL 0 OR NOT repeated STREQUAL "${single}${single}")
    message(FATAL_ERROR "streaming changed a field response or failed to stop on EOF")
endif()
