cmake_minimum_required(VERSION 2.8)

project(HF)

set(CMAKE_CXX_STANDARD 14)

# set up gsl
find_package(
        GSL 2.5
        REQUIRED
)

if( ${GSL_FOUND} )
    message("-- GSL library found")
    message("     GSL_VERSION: " ${GSL_VERSION_STRING})
    message("     GSL_LIBRARIES   : " ${GSL_LIBRARIES})
    message("     GSL_INCLUDE_DIRS: " ${GSL_INCLUDE_DIRS})
endif()

# directories to be searched for includes
#
include_directories( include )
include_directories( ${CMAKE_BINARY_DIR}/include )
include_directories(${GSL_INCLUDE_DIRS})

# Set up output directories
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# set up library
add_library(libHF
        src/basis.cpp
        src/gslextra.cpp
        src/gslprint.cpp
        src/integral.cpp
        src/RHF.cpp
        )


target_link_libraries(libHF GSL::gsl GSL::gslcblas)

add_executable(HF src/main.cpp)

target_link_libraries(HF libHF)