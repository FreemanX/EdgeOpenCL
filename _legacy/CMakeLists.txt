cmake_minimum_required(VERSION 3.12)

# External libraries
set(OpenCL_INCLUDE_DIR "PATH to OpenCL Library Include dir")
# Can be pulled from the smartphone if vendor provides OpenCL library for their GPU
set(OpenCL_LIBRARY "PATH to libOpenCL.so or libOpenCL.so") 
include_directories(${OpenCL_INCLUDE_DIR})
link_libraries(${OpenCL_LIBRARY})

# Compiling setttings
aux_source_directory(src DIR_LIB_SRC)
aux_source_directory(experiments DIR_EXAMPLE_SRC) 

include_directories(include)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/experiments)

add_library(heteroCompLib STATIC ${DIR_LIB_SRC})

file(GLOB APP_SOURCES experiments/*.cpp)
foreach(exampleSourceFile ${APP_SOURCES})
    string(REPLACE ".cpp" "" exampleName ${exampleSourceFile})
    string(REPLACE "/" ";" exampleName ${exampleName})
    list(REVERSE exampleName)
    list(GET exampleName 0 exampleName)
    message(STATUS "Example name: " ${exampleName})
    add_executable(${exampleName} ${exampleSourceFile})
    target_link_libraries(${exampleName} heteroCompLib)
    #target_link_libraries(${exampleName} external_libs)
endforeach(exampleSourceFile ${APP_SOURCES})


file(GLOB APP_SOURCES parboil/*.cpp)
foreach(exampleSourceFile ${APP_SOURCES})
    string(REPLACE ".cpp" "" exampleName ${exampleSourceFile})
    string(REPLACE "/" ";" exampleName ${exampleName})
    list(REVERSE exampleName)
    list(GET exampleName 0 exampleName)
    message(STATUS "Parboil app: " ${exampleName})
    add_executable(${exampleName} ${exampleSourceFile})
    target_link_libraries(${exampleName} heteroCompLib)
    #target_link_libraries(${exampleName} external_libs)
endforeach(exampleSourceFile ${APP_SOURCES})
