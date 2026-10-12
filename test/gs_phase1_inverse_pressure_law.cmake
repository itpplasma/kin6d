# CLI contract: pprime_law_file replaces the scalar pprime in the q_file mode.
file(MAKE_DIRECTORY "${WORK}")
execute_process(COMMAND "${DECK}" 2 --write-input WORKING_DIRECTORY "${WORK}"
    RESULT_VARIABLE status)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "inverse input generation failed")
endif()
# Normalized law, one cubic interval: pp(s) = P0*(1 + 0.2 s), ff unused.
set(P0 -1.4e5)
file(WRITE "${WORK}/pp.law"
"2 1
0.0 1.0 -1.0
0.0 1.0
${P0} -2.8e4 0.0 0.0 0.0 0.0 0.0 0.0
")
file(READ "${WORK}/input.nml" nml)
string(REGEX REPLACE "pprime= *[-+0-9.eE]+" "pprime=0.0, pprime_law_file='pp.law'"
    nml "${nml}")
file(WRITE "${WORK}/input.nml" "${nml}")
execute_process(COMMAND "${APP}" --config input.nml WORKING_DIRECTORY "${WORK}"
    RESULT_VARIABLE status OUTPUT_VARIABLE out ERROR_VARIABLE err TIMEOUT 120)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "variable-pressure inverse failed: ${out}${err}")
endif()
file(STRINGS "${WORK}/profile-law.dat" lines)
list(GET lines 3 row)
string(REGEX MATCHALL "[-+]?[0-9.]+[eE][-+][0-9]+" nums "${row}")
list(GET nums 0 pp0)
list(GET nums 1 pp1)
if(NOT pp0 MATCHES "^-1\\.4")
    message(FATAL_ERROR "recovered law lost the prescribed p' (${pp0})")
endif()
if(NOT pp1 MATCHES "^-2\\.8")
    message(FATAL_ERROR "recovered law lost the p' slope (${pp1})")
endif()
# A scalar pprime together with a pressure law is ambiguous and must fail.
file(READ "${WORK}/input.nml" nml)
string(REPLACE "pprime=0.0," "pprime=-1.0," nml "${nml}")
file(WRITE "${WORK}/conflict.nml" "${nml}")
execute_process(COMMAND "${APP}" --config conflict.nml WORKING_DIRECTORY "${WORK}"
    RESULT_VARIABLE status OUTPUT_QUIET ERROR_QUIET TIMEOUT 120)
if(status EQUAL 0)
    message(FATAL_ERROR "pprime with pprime_law_file accepted")
endif()
