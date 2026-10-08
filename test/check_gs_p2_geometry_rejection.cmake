execute_process(COMMAND "${GS_REFERENCE_ERROR}" "${MESH}" "${PARAMETERS}" affine-p2
    RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE diagnostic)
if(status EQUAL 0)
    message(FATAL_ERROR "Unsupported P2 midpoint geometry was accepted")
endif()
string(FIND "${diagnostic}" "P2 midpoint metadata is not native affine geometry" matched)
if(matched EQUAL -1)
    message(FATAL_ERROR "Wrong rejection diagnostic: ${diagnostic}")
endif()
