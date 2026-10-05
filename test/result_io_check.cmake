# F5 contract check, run by CTest with -P. Inputs: PROBE (test_result_io),
# CLASSICAL (classical_pair_demo), WORK (scratch directory), FORTIO_PIN.
# Oracles: FortIO's HDF5 reader inside the probe; system ncdump/h5dump when
# installed (NetCDF-C/HDF5-C are independent implementations).
file(REMOVE_RECURSE "${WORK}")
file(MAKE_DIRECTORY "${WORK}")

function(run_expect expect_ok pattern)
    execute_process(COMMAND ${ARGN} WORKING_DIRECTORY "${WORK}"
        RESULT_VARIABLE code OUTPUT_VARIABLE out ERROR_VARIABLE err)
    set(all "${out}${err}")
    if(expect_ok AND NOT code EQUAL 0)
        message(FATAL_ERROR "expected success: ${ARGN}\n${all}")
    endif()
    if(NOT expect_ok AND code EQUAL 0)
        message(FATAL_ERROR "expected rejection: ${ARGN}\n${all}")
    endif()
    if(NOT all MATCHES "${pattern}")
        message(FATAL_ERROR "output of ${ARGN} lacks /${pattern}/:\n${all}")
    endif()
    set(last_output "${all}" PARENT_SCOPE)
endfunction()

# Literal input with a comment, a partial array override and a string.
file(WRITE "${WORK}/literal.nml"
    "! probe override\n&probe scale = 2.5, weights(2) = 7, label = 'x y' /\n")
run_expect(TRUE "round trip PASS" "${PROBE}" --config literal.nml tag-a)
file(READ "${WORK}/output/test_result_io.resolved.nml" resolved_first)
if(NOT resolved_first MATCHES "SCALE= *2\\.5" OR
        NOT resolved_first MATCHES "WEIGHTS= *1\\.0+ *, *7\\.0+ *, *4\\.0+")
    message(FATAL_ERROR "resolved text lost overrides:\n${resolved_first}")
endif()
if(NOT resolved_first MATCHES "COUNT=3")
    message(FATAL_ERROR "resolved text lacks the default count:\n${resolved_first}")
endif()

find_program(NCDUMP ncdump)
if(NCDUMP)
    run_expect(TRUE "weights = 3, 21, 12 ;" "${NCDUMP}" -v weights,scale
        output/test_result_io.nc)
    foreach(pattern "scale = 2\\.5 ;" "double cube\\(cube_dim3, cube_dim2, cube_dim1\\)"
            "double grid\\(pair, w\\)" "scale:kin6d_status = \"certified-discrete\""
            "weights:units = \"m\"" ":kin6d_config_literal = \"! probe override"
            ":kin6d_config_resolved = \"&PROBE" ":fortio_revision = \"${FORTIO_PIN}\""
            ":kin6d_git_commit = \"[0-9a-f]+\"" ":created_utc = \"[0-9-]+T[0-9:]+Z\""
            ":compiler_version = \"" ":probe_argument = \"tag-a\"")
        if(NOT last_output MATCHES "${pattern}")
            message(FATAL_ERROR "ncdump lacks /${pattern}/:\n${last_output}")
        endif()
    endforeach()
    run_expect(TRUE "cube =" "${NCDUMP}" -v cube output/test_result_io.nc)
    # cube(i,j,k) = 2.5*(i + 10 j + 100 k); C order lists i fastest.
    string(REGEX REPLACE "[ \n]+" " " flat "${last_output}")
    if(NOT flat MATCHES "cube = 277\\.5, 280, 302\\.5, 305, 327\\.5, 330, 527\\.5")
        message(FATAL_ERROR "ncdump cube values:\n${last_output}")
    endif()
else()
    message(STATUS "ncdump not found; FortIO reader oracle only")
endif()
find_program(H5DUMP h5dump)
if(H5DUMP)
    run_expect(TRUE "\\(0\\): 3, 21, 12" "${H5DUMP}" -d /weights
        output/test_result_io.nc)
endif()

# Resolved text is itself a valid literal input and resolves to itself.
file(WRITE "${WORK}/resolved.nml" "${resolved_first}")
run_expect(TRUE "round trip PASS" "${PROBE}" --config=resolved.nml)
file(READ "${WORK}/output/test_result_io.resolved.nml" resolved_second)
if(NOT resolved_second STREQUAL resolved_first)
    message(FATAL_ERROR "resolve(resolve(x)) differs:\n${resolved_second}")
endif()

# Rejections: unknown entry, unknown group, malformed value, missing file.
file(WRITE "${WORK}/entry.nml" "&probe scal = 2.5 /\n")
run_expect(FALSE "invalid namelist &probe.*scal" "${PROBE}" --config entry.nml)
file(WRITE "${WORK}/group.nml" "&probe scale = 2.5 /\n&prob count = 2 /\n")
run_expect(FALSE "unknown namelist group &prob" "${PROBE}" --config group.nml)
file(WRITE "${WORK}/value.nml" "&probe count = 'two' /\n")
run_expect(FALSE "invalid namelist &probe" "${PROBE}" --config value.nml)
run_expect(FALSE "cannot open configuration" "${PROBE}" --config absent.nml)

# The classical pair rejects a shared &periodic_pair it cannot realize
# before any expensive work.
file(WRITE "${WORK}/pair.nml" "&periodic_pair coupling = 3.0 /\n")
run_expect(FALSE "periodic_pair mass/coupling differ" "${CLASSICAL}" --config pair.nml)
file(WRITE "${WORK}/legacy.nml" "&classical_pair length = 3.0 /\n")
run_expect(FALSE "invalid namelist &classical_pair" "${CLASSICAL}"
    --config legacy.nml)
# The box/time grid moved to the shared group; the old per-demo spelling
# is rejected instead of being silently ignored.
file(WRITE "${WORK}/quantum.nml" "&quantum_reference length = 3.0 /\n")
run_expect(FALSE "invalid namelist &quantum_reference.*length" "${QUANTUM}"
    --config quantum.nml)
