cmake_minimum_required(VERSION 3.16)
project(PRSSparseSolver LANGUAGES CXX)

if(POLICY CMP0144)
    cmake_policy(SET CMP0144 NEW)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

set(SOLVER_SOURCES
    src/matrix.cpp
    src/md.cpp
    src/permute.cpp
    src/symbolic.cpp
    src/numeric.cpp)

add_library(prs_solver STATIC ${SOLVER_SOURCES})
target_include_directories(prs_solver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(program src/main.cpp src/validate.cpp)
target_link_libraries(program PRIVATE prs_solver)

add_executable(bench
    benchmark/bench_superlu.cpp
    benchmark/superlu_adapter.cpp)
target_include_directories(bench PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/benchmark)
target_link_libraries(bench PRIVATE prs_solver)

if(MSVC)
    target_compile_options(program PRIVATE /utf-8)
    target_compile_options(bench PRIVATE /utf-8)
endif()

set(SUPERLU_ROOT "" CACHE PATH "SuperLU installation prefix")
set(SUPERLU_INCLUDE_DIR "" CACHE PATH "Directory containing slu_zdefs.h")
set(SUPERLU_LIBRARY "" CACHE FILEPATH "SuperLU library")
set(SUPERLU_BLAS_LIBRARIES "" CACHE STRING "BLAS libraries required by SuperLU")

if(SUPERLU_ROOT)
    list(PREPEND CMAKE_PREFIX_PATH "${SUPERLU_ROOT}")
endif()

find_package(superlu CONFIG QUIET)
if(NOT superlu_FOUND)
    find_package(SuperLU CONFIG QUIET)
endif()

if(TARGET superlu::superlu)
    set(SUPERLU_TARGET superlu::superlu)
elseif(TARGET SuperLU::SuperLU)
    set(SUPERLU_TARGET SuperLU::SuperLU)
elseif(TARGET superlu)
    set(SUPERLU_TARGET superlu)
endif()

if(NOT SUPERLU_TARGET)
    if(NOT SUPERLU_INCLUDE_DIR)
        find_path(SUPERLU_INCLUDE_DIR NAMES slu_zdefs.h
            HINTS
                "${SUPERLU_ROOT}/include"
                "${SUPERLU_ROOT}/include/superlu"
            PATHS
                /usr/include
                /usr/include/superlu
                /usr/local/include
                /usr/local/include/superlu)
    endif()

    if(NOT SUPERLU_LIBRARY)
        find_library(SUPERLU_LIBRARY NAMES superlu
            HINTS
                "${SUPERLU_ROOT}/lib"
                "${SUPERLU_ROOT}/lib64"
                "${SUPERLU_ROOT}/lib/x64"
            PATHS
                /usr/lib
                /usr/local/lib)
    endif()

    if(NOT SUPERLU_INCLUDE_DIR OR NOT SUPERLU_LIBRARY)
        message(FATAL_ERROR
            "SuperLU was not found. Set SUPERLU_ROOT, or set SUPERLU_INCLUDE_DIR and SUPERLU_LIBRARY.")
    endif()

    add_library(prs_superlu UNKNOWN IMPORTED)
    set_target_properties(prs_superlu PROPERTIES
        IMPORTED_LOCATION "${SUPERLU_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${SUPERLU_INCLUDE_DIR}")
    set(SUPERLU_TARGET prs_superlu)

    if(SUPERLU_BLAS_LIBRARIES)
        set_property(TARGET prs_superlu APPEND PROPERTY
            INTERFACE_LINK_LIBRARIES ${SUPERLU_BLAS_LIBRARIES})
    else()
        find_package(BLAS QUIET)
        if(BLAS_FOUND)
            set_property(TARGET prs_superlu APPEND PROPERTY
                INTERFACE_LINK_LIBRARIES BLAS::BLAS)
        else()
            find_library(PRS_BLAS_LIBRARY NAMES openblas blas
                HINTS
                    "${SUPERLU_ROOT}/lib"
                    "${SUPERLU_ROOT}/lib64")
            if(PRS_BLAS_LIBRARY)
                set_property(TARGET prs_superlu APPEND PROPERTY
                    INTERFACE_LINK_LIBRARIES "${PRS_BLAS_LIBRARY}")
            endif()
        endif()
    endif()
endif()

target_link_libraries(bench PRIVATE ${SUPERLU_TARGET})
