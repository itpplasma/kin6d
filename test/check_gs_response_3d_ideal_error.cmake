# Exercise the final-state report with eliminated resonant sheets and the
# field-aligned boundary sheet, independently of the zero-sheet unit oracle.
file(REMOVE_RECURSE "${WORK}")
file(MAKE_DIRECTORY "${WORK}")
file(WRITE "${WORK}/coils.txt"
    "1\nloop 4\n8 -1 -0.5\n8 1 -0.5\n8 1 0.5\n8 -1 0.5\n")
file(WRITE "${WORK}/currents.txt" "1000.0\n")
file(WRITE "${WORK}/input.nml"
    "&l4 eqdsk='${EQDSK}', coils='coils.txt', currents(1)='currents.txt', prefix='t', n=1, rat_num(1)=2, rat_den(1)=1, nrad=12, nth_max=64, order=2, nphi=32, nfine=256, np_sample=64, model='ideal', gmres_tol=1e-10, gmres_maxit=40 /\n")
execute_process(COMMAND "${APP}" input.nml WORKING_DIRECTORY "${WORK}"
    RESULT_VARIABLE status OUTPUT_VARIABLE out ERROR_VARIABLE err TIMEOUT 60)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "ideal response failed: ${out}${err}")
endif()
if(NOT out MATCHES "ideal response: info 0, GMRES iterations [1-9]")
    message(FATAL_ERROR "ideal GMRES did not converge: ${out}")
endif()
file(READ "${WORK}/t_res.dat" res)
if(NOT res MATCHES "# ideal_error phasing 1 status 0 frozen_defect_energy")
    message(FATAL_ERROR "missing successful final-state diagnostic: ${res}")
endif()
if(NOT res MATCHES "discretization unavailable-current-reconstruction sampling unavailable-current-sampling continuous unavailable-ideal-stability")
    message(FATAL_ERROR "missing separate unavailable error contributions: ${res}")
endif()
if(NOT res MATCHES "contraction directional-probe-not-upper-bound")
    message(FATAL_ERROR "directional contraction ratio mislabeled: ${res}")
endif()
