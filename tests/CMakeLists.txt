find_package(Catch2 3 REQUIRED)

add_executable(Tests
  GridBuilder.h
  GridBuilder.cpp
  TestGridBuilder.cpp
  TestFVM.cpp
)

target_link_libraries(Tests PRIVATE fluid_dynamics)
target_link_libraries(Tests PRIVATE Catch2::Catch2WithMain)

include(Catch)
catch_discover_tests(Tests)
