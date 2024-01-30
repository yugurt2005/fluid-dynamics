find_package(Catch2 3 REQUIRED)

add_executable(Tests
    FluidDynamicsTest.cpp
    FVMTest.cpp
    GridTest.cpp
    FaceTest.cpp
    SIMPLETest.cpp
    GridHelper.h
    GridHelper.cpp
)

target_link_libraries(Tests PRIVATE fluid_dynamics)
target_link_libraries(Tests PRIVATE Catch2::Catch2WithMain)

include(Catch)
catch_discover_tests(Tests)
