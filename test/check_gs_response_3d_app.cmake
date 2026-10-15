# Contract of the L4 response app: the discrete Ampere residual-energy
# indicator of its sheet solves and its claim status reach stdout and
# prefix_res.dat, with the continuous claim explicitly unavailable.
file(REMOVE_RECURSE "${WORK}")
file(MAKE_DIRECTORY "${WORK}")
# One horizontal loop; the coordinates only need to be a valid closed polyline.
set(coil "1\nloop 4\n8 0 0.5\n0 8 0.5\n-8 0 0.5\n0 -8 0.5\n")
file(WRITE "${WORK}/coils.txt" "${coil}")
file(WRITE "${WORK}/currents.txt" "1000.0\n")
file(WRITE "${WORK}/input.nml"
    "&l4 eqdsk='${EQDSK}', coils='coils.txt', currents(1)='currents.txt', prefix='t', n=1, rat_num(1)=2, rat_den(1)=1, nrad=12, nth_max=64, order=2, nphi=32, nfine=256, model='shield' /\n")
execute_process(COMMAND "${APP}" input.nml WORKING_DIRECTORY "${WORK}"
    RESULT_VARIABLE status OUTPUT_VARIABLE out ERROR_VARIABLE err TIMEOUT 120)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "gs_response_3d failed: ${out}${err}")
endif()
if(NOT out MATCHES "residual-energy indicator [^\n]*solves 1, without indicator 0; claim indicator, continuous error unavailable-exterior-remainder")
    message(FATAL_ERROR "missing Ampere indicator status in stdout: ${out}")
endif()
file(READ "${WORK}/t_res.dat" res)
if(NOT res MATCHES "# ampere_discrete_energy_error_max +([0-9.E+-]+) sheet_solves 1 without_indicator 0 claim indicator")
    message(FATAL_ERROR "missing Ampere indicator in t_res.dat: ${res}")
endif()
set(indicator "${CMAKE_MATCH_1}")
if(NOT indicator LESS 1e-6)
    message(FATAL_ERROR "algebraic Ampere residual energy ${indicator} not at solver level")
endif()
if(NOT res MATCHES "continuous_claim unavailable-exterior-remainder")
    message(FATAL_ERROR "continuous claim not marked unavailable")
endif()
