execute_process(COMMAND "${GS_REFERENCE_ERROR}" "${MESH}" "${PARAMETERS}"
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE diagnostic)
if(result EQUAL 0)
    message(FATAL_ERROR "P1 reference checker accepted a quadratic mesh")
endif()
if(NOT diagnostic MATCHES "reference enclosure supports P1 mesh-v1 only")
    message(FATAL_ERROR "Unexpected rejection: ${diagnostic}")
endif()
