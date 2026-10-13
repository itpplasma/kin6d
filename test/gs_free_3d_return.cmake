# NL-3a on the gs_free fixture: the 3D free boundary at N = 0 must reproduce
# the 2D free-boundary solution (gs_free, P3, refit) and return to it from
# offset starts (+-15 mm in R and Z, +2 % elongation). Gauge-invariant
# measures: largest distance of the mapped Gamma_p nodes from the 2D boundary
# (psi/|grad psi| of the 2D result) and the interior surface defect.
file(REMOVE_RECURSE "${WORK}")
file(MAKE_DIRECTORY "${WORK}/ref")
foreach(name boundary.dat coils.txt profile.dat limiter.dat)
    file(COPY "${FIXTURES}/${name}" DESTINATION "${WORK}/ref")
endforeach()
file(WRITE "${WORK}/ref/input.nml" "&gs_free boundary_file='boundary.dat', boundary_nodes=96, degree=3, coil_file='coils.txt', profile_file='profile.dat', plasma_ellipse=3.0,0.7,1.0, chain_nodes=64, lcfs_mode='limiter', points_file='limiter.dat', axis_guess=3.0,0.0, search_radius=0.5, seed_amplitude=0.5, max_refits=6, currents=-10000.0,0.0, newton_tolerance=1e-11 /\n")
execute_process(COMMAND "${GS_FREE}" --config input.nml WORKING_DIRECTORY "${WORK}/ref"
    RESULT_VARIABLE status OUTPUT_VARIABLE out ERROR_VARIABLE err)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "gs_free reference failed: ${out}${err}")
endif()

find_package(Python3 REQUIRED COMPONENTS Interpreter)
function(summary_value file name result)
    file(STRINGS "${file}" lines REGEX "^${name} ")
    string(REGEX REPLACE "^${name} +" "" value "${lines}")
    set(${result} "${value}" PARENT_SCOPE)
endfunction()
function(require_below what value limit)
    execute_process(COMMAND "${Python3_EXECUTABLE}" -I -c
        "import sys; sys.exit(0 if abs(float(sys.argv[1])) <= float(sys.argv[2]) else 1)"
        "${value}" "${limit}" RESULT_VARIABLE below)
    if(NOT below EQUAL 0)
        message(FATAL_ERROR "${what}: ${value} exceeds ${limit}")
    endif()
endfunction()

set(common "mesh_file='../ref/mesh.dat', law_file='../ref/profile-law.dat', summary_file='../ref/summary.dat', coil_file='../ref/coils.txt', boundary_file='../ref/boundary.dat', boundary_nodes=64, degree=3, chain_nodes=48, nfp=1, ntor=0, max_iterations=150, verbose=.false.")
foreach(case "r0;offset_r=0.0" "rp;offset_r=0.015" "rm;offset_r=-0.015" "zp;offset_z=0.015" "zm;offset_z=-0.015" "el;elongation=0.02"
    "gr0;offset_r=0.0, tangential_gauge=100.0"
    "gzm;offset_z=-0.015, tangential_gauge=100.0, rigid_vertical=.true.")
    list(GET case 0 name)
    list(GET case 1 offset)
    file(MAKE_DIRECTORY "${WORK}/${name}")
    file(WRITE "${WORK}/${name}/in.nml" "&gs_free_3d ${common}, ${offset} /\n")
    execute_process(COMMAND "${APP}" --config in.nml WORKING_DIRECTORY "${WORK}/${name}"
        RESULT_VARIABLE status OUTPUT_VARIABLE out ERROR_VARIABLE err)
    if(NOT status EQUAL 0)
        message(FATAL_ERROR "gs_free_3d ${name} failed: ${out}${err}")
    endif()
    message(STATUS "${name}: ${out}")
    set(summary "${WORK}/${name}/summary-3d.dat")
    summary_value("${summary}" boundary_distance_max_m distance)
    summary_value("${summary}" surface_defect_max_m defect)
    require_below("${name} boundary distance (m)" "${distance}" "${DISTANCE_LIMIT}")
    require_below("${name} surface defect (m)" "${defect}" "${DEFECT_LIMIT}")
endforeach()
