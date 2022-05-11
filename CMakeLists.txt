cmake_minimum_required(VERSION 3.12)
# Project settings
project(EDCL) #Edge OpenCL
set (CMAKE_CXX_STANDARD 11)

# General cmake settings
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/examples)


# Android related settings
set(CMAKE_ANDROID_STANDALONE_TOOLCHAIN "$ENV{HOME}/toolchains/ndk_r17c")
include_directories(${CMAKE_ANDROID_STANDALONE_TOOLCHAIN}/lib64/clang/6.0.2/include/ ${CMAKE_ANDROID_STANDALONE_TOOLCHAIN}/sysroot/usr/include)
link_directories(${CMAKE_ANDROID_STANDALONE_TOOLCHAIN}/sysroot/usr/lib)

# EDCL Library
include_directories(include)
aux_source_directory(src LIBEDCL)
add_library(EDCL STATIC ${LIBEDCL})

# Include OpenCL Kernels
include_directories(kernels)

# Tests
aux_source_directory(tests TESTS)

file(GLOB APP_SOURCES tests/*.cpp)
foreach(exampleSourceFile ${APP_SOURCES})
    string(REPLACE ".cpp" "" exampleName ${exampleSourceFile})
    string(REPLACE "/" ";" exampleName ${exampleName})
    list(REVERSE exampleName)
    list(GET exampleName 0 exampleName)
    message(STATUS "Tests: " ${exampleName})
    add_executable(${exampleName} ${exampleSourceFile})
    target_link_libraries(${exampleName} EDCL)
endforeach(exampleSourceFile ${APP_SOURCES})

# Polybench-gpu 
aux_source_directory(polybench-gpu POLY)
file(GLOB APP_SOURCES polybench-gpu/*.cpp)
foreach(exampleSourceFile ${APP_SOURCES})
    string(REPLACE ".cpp" "" exampleName ${exampleSourceFile})
    string(REPLACE "/" ";" exampleName ${exampleName})
    list(REVERSE exampleName)
    list(GET exampleName 0 exampleName)
    message(STATUS "Polybench: " ${exampleName})
    add_executable(${exampleName} ${exampleSourceFile})
    target_link_libraries(${exampleName} EDCL)
endforeach(exampleSourceFile ${APP_SOURCES})
