#include <algorithm>

#include <catch2/catch_test_macros.hpp>

#include "Eigen/Core"
#include "Eigen/Dense"

#include "../src/FVM.h"
#include "../src/Grid.h"

using Eigen::VectorXd;

TEST_CASE("Grid")
{
  int n = 4;

  // 0 1
  // 3 2
  vector<std::pair<int, int>> connections = {
      {0, 1},
      {1, 2},
      {2, 3},
      {3, 0},
  };
  vector<double> cellX = {0, 1, 1, 0};
  vector<double> cellY = {0, 0, 1, 1};
  vector<std::pair<int, Vector2d>> walls = {};

  Grid grid(n, connections, cellX, cellY, walls);

  SECTION("Grid getAdjacents: Basic")
  {
    REQUIRE(grid.getAdjacents(0) == vector<int>{1, 3});
    REQUIRE(grid.getAdjacents(1) == vector<int>{0, 2});
    REQUIRE(grid.getAdjacents(2) == vector<int>{1, 3});
    REQUIRE(grid.getAdjacents(3) == vector<int>{2, 0});
  }

  SECTION("Grid getAdjDx: Basic") {
    VectorXd dx0 = grid.getAdjDx(0);
    VectorXd dx1 = grid.getAdjDx(1);
    VectorXd dx2 = grid.getAdjDx(2);
    VectorXd dx3 = grid.getAdjDx(3);

    REQUIRE(dx0.rows() == 2);

    REQUIRE(dx0(0) == +1);
    REQUIRE(dx0(1) == 0);

    REQUIRE(dx1(0) == -1);
    REQUIRE(dx1(1) == 0);

    REQUIRE(dx2(0) == 0);
    REQUIRE(dx2(1) == -1);

    REQUIRE(dx3(0) == +1);
    REQUIRE(dx3(1) == 0);
  }
}