#include <catch2/catch_test_macros.hpp>

#include <iostream>

#include "Eigen/Core"
#include "Eigen/Dense"

#include "../src/SIMPLE.h"
#include "../src/FVM.h"
#include "../src/Grid.h"
#include "GridHelper.h"

TEST_CASE("SIMPLE loop: Happy Test")
{
  int n = 6;

  Parameters parameters;

  //          0 1 2
  //          - - -
  // 3 4 5    . . . | 1
  // 0 1 2    . . . | 0
  Grid grid = buildRectangularGrid(2, 3, 1);

  FVM fvm(grid);

  SIMPLE simple(parameters, fvm);

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

    State result = simple.loop(state, 1);
  }
}

TEST_CASE("SIMPLE loop: Convergence")
{
  int n = 6;

  Parameters parameters;

  //          0 1 2
  //          - - -
  // 3 4 5    . . . | 1
  // 0 1 2    . . . | 0
  Grid grid = buildRectangularGrid(2, 3, 1);

  FVM fvm(grid);

  SIMPLE simple(parameters, fvm);

  SECTION("SIMPLE loop: Solved")
  {
    VectorXd u(n);
    VectorXd v(n);
    VectorXd p(n);
    u << 1, 1, 1, 1, 1, 1;
    v << 0, 0, 0, 0, 0, 0;
    p << 1, 0, 0, 1, 0, 0;

    State state = {
        u,
        v,
        p,
        VectorXd::Zero(n),
        VectorXd::Zero(n),
    };

    State result = simple.loop(state, 10);
  }
}