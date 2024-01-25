#include <catch2/catch_test_macros.hpp>

#include <iostream>

#include "Eigen/Core"
#include "Eigen/Dense"

#include "../src/SIMPLE.h"
#include "../src/FVM.h"
#include "../src/Grid.h"

TEST_CASE("SIMPLE calculateVelocity: Laminar")
{
  int n = 6;

  Parameters parameters;
  //          0 1 3
  //          - - -
  // 3 4 5    . . . | 3
  // 0 1 2    . . . | 1
  vector<std::pair<int, int>> connections = {
      {0, 1},
      {1, 2},
      {3, 4},
      {4, 5},
      {0, 3},
      {1, 4},
      {2, 5},
  };
  vector<double> cellX = {0, 1, 3, 0, 1, 3};
  vector<double> cellY = {1, 1, 1, 3, 3, 3};
  vector<std::pair<int, Vector2d>> walls = {
      {0, Vector2d(1, 0)},
      {3, Vector2d(1, 0)},
  };

  Grid grid(n, connections, cellX, cellY, walls);

  FVM fvm(n, grid);

  SIMPLE simple(parameters, fvm, grid);

  VectorXd u(n);
  VectorXd v(n);
  VectorXd p(n);
  u << 1, 1, 1, 1, 1, 1;
  v << 0, 0, 0, 0, 0, 0;
  p << 2, 1, 0, 2, 1, 0;

  State state = {
      u,
      v,
      p,
      VectorXd::Zero(1),
      VectorXd::Zero(1),
  };

  simple.calculateVelocity(state);

  VectorXd &newU = simple.getNewU();
  VectorXd &newV = simple.getNewV();
  for (int i = 0; i < n; i++)
  {
    REQUIRE(newV(i) == 0);
    REQUIRE(newU(i) == 1);
  }
}

TEST_CASE("SIMPLE loop: Convergence")
{
  int n = 6;

  Parameters parameters;

  //          0 1 3
  //          - - -
  // 3 4 5    . . . | 3
  // 0 1 2    . . . | 1
  vector<std::pair<int, int>> connections = {
      {0, 1},
      {1, 2},
      {3, 4},
      {4, 5},
      {0, 3},
      {1, 4},
      {2, 5},
  };
  vector<double> cellX = {0, 1, 3, 0, 1, 3};
  vector<double> cellY = {1, 1, 1, 3, 3, 3};
  vector<std::pair<int, Vector2d>> walls = {
      {0, Vector2d(1, 0)},
      {3, Vector2d(1, 0)},
  };

  Grid grid(n, connections, cellX, cellY, walls);

  FVM fvm(n, grid);

  SIMPLE simple(parameters, fvm, grid);

  SECTION("SIMPLE loop: Solved")
  {
    VectorXd u(n);
    VectorXd v(n);
    u << 1, 1, 1, 1, 1, 1;
    v << 0, 0, 0, 0, 0, 0;

    State state = {
        u,
        v,
        VectorXd::Random(n),
        VectorXd::Zero(n),
        VectorXd::Zero(n),
    };

    State result = simple.loop(state, 6);
  }

  SECTION("SIMPLE loop: Basic")
  {
    VectorXd u(n);
    VectorXd v(n);
    u << 1, -1, -1, 1, -1, -1;
    v << 0, 0, 0, 0, 0, 0;

    State state = {
        u,
        v,
        VectorXd::Zero(n),
        VectorXd::Zero(n),
        VectorXd::Zero(n),
    };

    // TODO
  }

  SECTION("SIMPLE loop: Mixed")
  {
    VectorXd u(n);
    VectorXd v(n);
    u << 1, -1, -1, 1, -1, -1;
    v << 0, 1, -1, 0, -1, 1;

    State state = {
        u,
        v,
        VectorXd::Zero(n),
        VectorXd::Zero(n),
        VectorXd::Zero(n),
    };

    // TODO
  }
}