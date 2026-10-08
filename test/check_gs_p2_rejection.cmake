execute_process(COMMAND "${GS_REFERENCE_ERROR}" "${MESH}" "${PARAMETERS}"
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE diagnostic)
if(result EQUAL 0)
    message(FATAL_ERROR "P1 reference checker accepted a quadratic mesh")
endif()
if(NOT diagnostic MATCHES "reference enclosure supports P1 mesh-v1 only")
    message(FATAL_ERROR "Unexpected rejection: ${diagnostic}")
endif()

# A forged header must not let a P1 reader truncate six-node connectivity.
file(READ "${MESH}" mesh_text)
string(REGEX REPLACE "^[^\n]*" "# kin6d-gs-mesh-v1: forged P2 connectivity"
    forged_mesh_text "${mesh_text}")
set(forged_mesh "${MESH}.forged-v1")
file(WRITE "${forged_mesh}" "${forged_mesh_text}")
execute_process(COMMAND "${GS_REFERENCE_ERROR}" "${forged_mesh}" "${PARAMETERS}"
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE diagnostic)
if(result EQUAL 0)
    message(FATAL_ERROR "P1 reference checker truncated six-node connectivity")
endif()
if(NOT diagnostic MATCHES "extra connectivity entries")
    message(FATAL_ERROR "Unexpected forged-header rejection: ${diagnostic}")
endif()
