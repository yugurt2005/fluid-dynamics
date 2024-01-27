cmake_minimum_required(VERSION 3.27)
project(fluid_dynamics)

set(CMAKE_CXX_STANDARD 17)

add_library(fluid_dynamics STATIC
  src/FluidDynamics.h
  src/FluidDynamics.cpp
  src/SIMPLE.h
  src/SIMPLE.cpp
  src/Parameters.h
  src/State.h
  src/State.cpp
  src/FVM.h
  src/FVM.cpp
  src/Grid.h
  src/Grid.cpp
  src/Face.h
  src/Face.cpp
  interfaces/IGrid.h
  )

find_package(Eigen3 REQUIRED)
target_link_libraries(fluid_dynamics
  Eigen3::Eigen)

add_subdirectory(tests)
