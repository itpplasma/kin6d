# Dependency provider for the dev-fast preset (CMAKE_PROJECT_TOP_LEVEL_INCLUDES).
# FetchContent_MakeAvailable() sources are cloned once per name and pin into
# KIN6D_SHARED_DEPS_DIR and shared by all worktrees and build trees; every build
# tree keeps its own binary directory. A pin change gets a new directory, so a
# shared checkout is never rewritten under a running build. Dependencies that
# KIN6D populates without MakeAvailable (fortnum, fortfem) stay per build tree.
set(KIN6D_SHARED_DEPS_DIR "$ENV{HOME}/.cache/kin6d/deps" CACHE PATH
    "Shared FetchContent source checkouts, keyed by name and pin")

macro(kin6d_shared_deps_provide method name)
    cmake_parse_arguments(_kin6d_dep "" "GIT_REPOSITORY;GIT_TAG;SOURCE_SUBDIR;BINARY_DIR"
        "" ${ARGN})
    if(_kin6d_dep_GIT_REPOSITORY AND _kin6d_dep_GIT_TAG)
        string(TOLOWER "${name}" _kin6d_dep_lower)
        set(_kin6d_dep_src "${KIN6D_SHARED_DEPS_DIR}/${_kin6d_dep_lower}-${_kin6d_dep_GIT_TAG}")
        set(_kin6d_dep_args ${ARGN})
        foreach(_kin6d_dep_drop SOURCE_DIR BINARY_DIR SOURCE_SUBDIR)
            list(FIND _kin6d_dep_args ${_kin6d_dep_drop} _kin6d_dep_at)
            if(_kin6d_dep_at GREATER -1)
                list(REMOVE_AT _kin6d_dep_args ${_kin6d_dep_at})
                list(REMOVE_AT _kin6d_dep_args ${_kin6d_dep_at})
            endif()
        endforeach()
        list(REMOVE_ITEM _kin6d_dep_args EXCLUDE_FROM_ALL SYSTEM)
        file(MAKE_DIRECTORY "${KIN6D_SHARED_DEPS_DIR}")
        file(LOCK "${_kin6d_dep_src}.lock" GUARD PROCESS TIMEOUT 1800)
        if(NOT EXISTS "${_kin6d_dep_src}/.kin6d-complete")
            file(REMOVE_RECURSE "${_kin6d_dep_src}")
            FetchContent_Populate(kin6d_shared_${_kin6d_dep_lower} QUIET
                SOURCE_DIR "${_kin6d_dep_src}"
                SUBBUILD_DIR "${CMAKE_BINARY_DIR}/_deps/${_kin6d_dep_lower}-shared-subbuild"
                BINARY_DIR "${CMAKE_BINARY_DIR}/_deps/${_kin6d_dep_lower}-shared-unused"
                ${_kin6d_dep_args})
            file(TOUCH "${_kin6d_dep_src}/.kin6d-complete")
        endif()
        file(LOCK "${_kin6d_dep_src}.lock" RELEASE)
        FetchContent_SetPopulated(${name} SOURCE_DIR "${_kin6d_dep_src}"
            BINARY_DIR "${_kin6d_dep_BINARY_DIR}")
        if(_kin6d_dep_SOURCE_SUBDIR)
            add_subdirectory("${_kin6d_dep_src}/${_kin6d_dep_SOURCE_SUBDIR}"
                "${_kin6d_dep_BINARY_DIR}")
        else()
            add_subdirectory("${_kin6d_dep_src}" "${_kin6d_dep_BINARY_DIR}")
        endif()
    endif()
endmacro()

cmake_language(SET_DEPENDENCY_PROVIDER kin6d_shared_deps_provide
    SUPPORTED_METHODS FETCHCONTENT_MAKEAVAILABLE_SERIAL)
