cmake_minimum_required(VERSION 3.16.3)

project(openFDTD
  DESCRIPTION "A free and open source FDTD solver"
  LANGUAGES C CXX)

set(CMAKE_CXX_FLAGS "-Wall -Wpedantic -Wextra -fopenmp -lglfw3 -lGL -lX11 -lpthread -lXrandr -lXi -ldl")

set(CMAKE_THREAD_LIBS_INIT "-lpthread")
set(CMAKE_HAVE_THREADS_LIBRARY 1)
set(CMAKE_USE_WIN32_THREADS_INIT 0)
set(CMAKE_USE_PTHREADS_INIT 1)
set(THREADS_PREFER_PTHREAD_FLAG ON)


set(CMAKE_SOURCE_DIR ..)
include_directories(${CMAKE_SOURCE_DIR}/include)
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

find_package(glfw3 REQUIRED)
find_package(OpenGL REQUIRED)

file(GLOB SOURCES "src/*.cpp")

add_executable(openFDTD ${SOURCES})
target_link_libraries(openFDTD glfw)
target_link_libraries(openFDTD OpenGL::GL)
target_link_libraries(openFDTD ${CMAKE_DL_LIBS})
