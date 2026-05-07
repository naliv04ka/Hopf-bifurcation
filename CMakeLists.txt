cmake_minimum_required(VERSION 3.10)
project(CavityFlow)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -Wall -Wextra")

# Source files
set(SOURCES
    src/main.cpp
    src/grid.cpp
    src/field.cpp
    src/pimple_solver.cpp
    src/fft_analyzer.cpp
    src/post_processor.cpp
)

# Include directory
include_directories(src)

# Executable
add_executable(cavity_flow ${SOURCES})

# Math library
target_link_libraries(cavity_flow m)

# Filesystem library (required for GCC < 9)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
    target_link_libraries(cavity_flow stdc++fs)
endif()
