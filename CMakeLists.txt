cmake_minimum_required(VERSION 3.27)
project(fluid_dynamics)

set(CMAKE_CXX_STANDARD 17)

add_library(fluid_dynamics STATIC
  src/Mesh.h
  src/Mesh.cpp
  src/Cell.h
  src/Cell.cpp
  src/Edge.h
  src/Edge.cpp
  src/Debug.h
  src/Debug.cpp
  )

find_package(Eigen3 REQUIRED)
target_link_libraries(fluid_dynamics
  Eigen3::Eigen)

add_subdirectory(tests)