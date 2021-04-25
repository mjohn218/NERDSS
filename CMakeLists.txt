cmake_minimum_required(VERSION 3.0)

project(nerdss)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -fsanitize=address -fno-omit-frame-pointer")

set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(default_build_type "Release")
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to default, '${default_build_type}'.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build." FORCE)
endif()

if(${PROJECT_SOURCE_DIR} STREQUAL ${PROJECT_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed. Please use a different directory (e.g. 'build').")
endif()

# Get system information
cmake_host_system_information(RESULT NUM_CPUS QUERY NUMBER_OF_PHYSICAL_CORES)
cmake_host_system_information(RESULT NUM_THREADS QUERY NUMBER_OF_LOGICAL_CORES)

# Print system information
message(STATUS "We are on a ${CMAKE_SYSTEM_NAME} system")
message(STATUS "The host processor is ${CMAKE_HOST_SYSTEM_PROCESSOR} with ${NUM_CPUS} physical cores and ${NUM_THREADS} logical cores.")

include_directories(include)

#add_definitions(-D_DEBUG) # comment this line if no need _DEBUG

file(GLOB SOURCES "src/*/*.cpp")
add_executable(nerdss ${SOURCES} EXEs/nerdss.cpp)

# Set up external libraries
find_package(GSL REQUIRED)
target_link_libraries(nerdss GSL::gsl GSL::gslcblas)

# Set up header directories
include_directories(include $(GSL_INCLUDE_DIR))