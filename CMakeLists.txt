cmake_minimum_required(VERSION 3.15)
project(AxiomR
    VERSION 0.1.0
    DESCRIPTION "CPU based renderer"
    LANGUAGES CXX
)

message(STATUS "Available configurations: ${CMAKE_CONFIGURATION_TYPES}")
include(FetchContent)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compilation flags
if(MSVC)
else()
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# Detect architecture (example)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|amd64")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64")
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "arm")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=armv7-a") 
endif()

# Output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Options
option(BUILD_TESTS "Build test programs" ON)

# Dependencies
add_subdirectory(external)

# Main executable
add_subdirectory(src)

# Tests
if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
