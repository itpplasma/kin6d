# Dependency provider for the dev-fast preset (CMAKE_PROJECT_TOP_LEVEL_INCLUDES).
# FetchContent_MakeAvailable() sources are cloned once per name and pin into
# MILCA_SHARED_DEPS_DIR and shared by all worktrees and build trees; every build
# tree keeps its own binary directory. A pin change gets a new directory, so a
# shared checkout is never rewritten under a running build. Dependencies that
# MILCA populates without MakeAvailable (fortnum, fortfem) stay per build tree.
set(MILCA_SHARED_DEPS_DIR "$ENV{HOME}/.cache/milca/deps" CACHE PATH
    "Shared FetchContent source checkouts, keyed by name and pin")

macro(milca_shared_deps_provide method name)
    cmake_parse_arguments(_milca_dep "" "GIT_REPOSITORY;GIT_TAG;SOURCE_SUBDIR;BINARY_DIR"
        "" ${ARGN})
    if(_milca_dep_GIT_REPOSITORY AND _milca_dep_GIT_TAG)
        string(TOLOWER "${name}" _milca_dep_lower)
        set(_milca_dep_src "${MILCA_SHARED_DEPS_DIR}/${_milca_dep_lower}-${_milca_dep_GIT_TAG}")
        set(_milca_dep_args ${ARGN})
        foreach(_milca_dep_drop SOURCE_DIR BINARY_DIR SOURCE_SUBDIR)
            list(FIND _milca_dep_args ${_milca_dep_drop} _milca_dep_at)
            if(_milca_dep_at GREATER -1)
                list(REMOVE_AT _milca_dep_args ${_milca_dep_at})
                list(REMOVE_AT _milca_dep_args ${_milca_dep_at})
            endif()
        endforeach()
        list(REMOVE_ITEM _milca_dep_args EXCLUDE_FROM_ALL SYSTEM)
        file(MAKE_DIRECTORY "${MILCA_SHARED_DEPS_DIR}")
        file(LOCK "${_milca_dep_src}.lock" GUARD PROCESS TIMEOUT 1800)
        if(NOT EXISTS "${_milca_dep_src}/.milca-complete")
            file(REMOVE_RECURSE "${_milca_dep_src}")
            FetchContent_Populate(milca_shared_${_milca_dep_lower} QUIET
                SOURCE_DIR "${_milca_dep_src}"
                SUBBUILD_DIR "${CMAKE_BINARY_DIR}/_deps/${_milca_dep_lower}-shared-subbuild"
                BINARY_DIR "${CMAKE_BINARY_DIR}/_deps/${_milca_dep_lower}-shared-unused"
                ${_milca_dep_args})
            file(TOUCH "${_milca_dep_src}/.milca-complete")
        endif()
        file(LOCK "${_milca_dep_src}.lock" RELEASE)
        if(NOT _milca_dep_BINARY_DIR)
            set(_milca_dep_BINARY_DIR "${CMAKE_BINARY_DIR}/_deps/${_milca_dep_lower}-build")
        endif()
        FetchContent_SetPopulated(${name} SOURCE_DIR "${_milca_dep_src}"
            BINARY_DIR "${_milca_dep_BINARY_DIR}")
        # Like FetchContent_MakeAvailable: a SOURCE_SUBDIR without CMakeLists.txt
        # only populates (used for libneo).
        set(_milca_dep_cmake_dir "${_milca_dep_src}")
        if(_milca_dep_SOURCE_SUBDIR)
            set(_milca_dep_cmake_dir "${_milca_dep_src}/${_milca_dep_SOURCE_SUBDIR}")
        endif()
        if(EXISTS "${_milca_dep_cmake_dir}/CMakeLists.txt")
            add_subdirectory("${_milca_dep_cmake_dir}" "${_milca_dep_BINARY_DIR}")
        endif()
    endif()
endmacro()

cmake_language(SET_DEPENDENCY_PROVIDER milca_shared_deps_provide
    SUPPORTED_METHODS FETCHCONTENT_MAKEAVAILABLE_SERIAL)
