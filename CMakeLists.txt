cmake_minimum_required(VERSION 3.10)
project(pickett CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Default build type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release)
endif()

# For multi-config generators (like Visual Studio), set configurations
if(CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_CONFIGURATION_TYPES Release Debug RelWithDebInfo MinSizeRel)
endif()

# Compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -Wall -DNDEBUG")
set(CMAKE_CXX_FLAGS_MINSIZEREL "-Os -Wall -DNDEBUG")

# Find BLAS (optional)
find_package(BLAS)
if(BLAS_FOUND)
    set(BLAS_LIBRARIES ${BLAS_LIBRARIES})
    message(STATUS "Found BLAS: ${BLAS_LIBRARIES}")
else()
    message(STATUS "BLAS not found, using supplied routines")
    set(BLAS_LIBRARIES "")
endif()

# On Windows, if using MSVC, ensure C++11 support
if(MSVC)
    add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++11")
endif()

# Include directories
include_directories(.)

# Library: splib
set(SPLIB_SOURCES
    ulib.cpp
    cnjj.cpp
    slibgcc.cpp
    catutil.cpp
    lsqfit.cpp
    dblas.cpp
)
add_library(splib STATIC ${SPLIB_SOURCES})

# Executables
add_executable(spfit calfit.cpp subfit.cpp spinv.cpp spinit.cpp)
target_link_libraries(spfit splib ${BLAS_LIBRARIES} m)

add_executable(spcat calcat.cpp sortsub.cpp spinv.cpp spinit.cpp)
target_link_libraries(spcat splib ${BLAS_LIBRARIES} m)

add_executable(calmrg calmrg.cpp)
target_link_libraries(calmrg splib ${BLAS_LIBRARIES} m)

add_executable(dpfit calfit.cpp subfit.cpp dpi.cpp)
target_link_libraries(dpfit splib ${BLAS_LIBRARIES} m)

add_executable(dpcat calcat.cpp sortsub.cpp dpi.cpp)
target_link_libraries(dpcat splib ${BLAS_LIBRARIES} m)

# Additional executables
add_executable(iamcalc iamcalc.cpp ftran.cpp readopt.cpp)
target_link_libraries(iamcalc splib ${BLAS_LIBRARIES} m)

add_executable(moiam moiam.cpp ftran.cpp)
target_link_libraries(moiam splib ${BLAS_LIBRARIES} m)

add_executable(stark stark.cpp)
target_link_libraries(stark splib ${BLAS_LIBRARIES} m)

add_executable(termval termval.cpp)
target_link_libraries(termval splib ${BLAS_LIBRARIES} m)

add_executable(sortn sortn.cpp sortsub.cpp)
target_link_libraries(sortn splib ${BLAS_LIBRARIES} m)

add_executable(calbak calbak.cpp)
target_link_libraries(calbak splib ${BLAS_LIBRARIES} m)

