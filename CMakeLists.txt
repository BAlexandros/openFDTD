cmake_minimum_required(VERSION 3.16.3)

project(openFDTD
  DESCRIPTION "A free and open source FDTD solver"
  LANGUAGES CXX)

set(CMAKE_CXX_FLAGS "-Wall -Wpedantic -Wextra -fopenmp")

set(CMAKE_SOURCE_DIR ..)
include_directories(${CMAKE_SOURCE_DIR}/include)
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

file(GLOB SOURCES "src/*.cpp")

add_executable(openFDTD ${SOURCES})
