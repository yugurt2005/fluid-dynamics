#include <catch2/catch_test_macros.hpp>

#include "Eigen/Core"
#include "Eigen/Dense"

#include "../src/FVM.h"
#include "../src/Grid.h"

using Eigen::VectorXd;

TEST_CASE("FVM calcLeastSquaresGradient: Happy Test", "[Gradient]")
{
  int n = 4;

  VectorXd dxs(n);
  VectorXd dys(n);

  dxs << 0, 1, 0, -1;
  dys << 1, 0, -1, 0;

  VectorXd values(n);
  values << 1, 1, -1, -1;

  VectorXd result = FVM::calcLeastSquaresGradient(n, dxs, dys) * values;
  double dx = result(0);
  double dy = result(1);

  REQUIRE(dx == 1);
  REQUIRE(dy == 1);
}

TEST_CASE("FVM calcLeastSquaresGradient: Inconsistent", "[Gradient]")
{
  int n = 4;

  VectorXd dxs(n);
  VectorXd dys(n);

  dxs << 0, 1, 0, -1;
  dys << 1, 0, -1, 0;

  VectorXd values(n);
  values << 1, 1, 1, 1;

  VectorXd result = FVM::calcLeastSquaresGradient(n, dxs, dys) * values;
  double dx = result(0);
  double dy = result(1);

  REQUIRE(dx == 0);
  REQUIRE(dy == 0);
}

TEST_CASE("FVM calcLeastSquaresGradient: Weights", "[Gradient]")
{
  int n = 4;

  VectorXd dxs(n);
  VectorXd dys(n);

  dxs << 0, 1, 0, -2;
  dys << 1, 0, -2, 0;

  VectorXd values(n);
  values << 1, 1, -1, -1;

  VectorXd result = FVM::calcLeastSquaresGradient(n, dxs, dys) * values;
  double dx = result(0);
  double dy = result(1);

  REQUIRE(dx == 0.75);
  REQUIRE(dy == 0.75);
}

TEST_CASE("FVM createGradientMatrix: Happy Test", "[Gradient]")
{
  int n = 4;

  // 0 1    1 2
  // 3 2    0 3
  vector<std::pair<int, int>> connections = {
      {0, 1},
      {1, 2},
      {2, 3},
      {3, 0},
  };
  vector<double> cellX = {0, 1, 1, 0};
  vector<double> cellY = {1, 1, 0, 0};
  vector<std::pair<int, Vector2d>> walls = {};

  Grid grid(n, connections, cellX, cellY, walls);

  FVM fvm(n, grid);

  const SparseMatrix<double> &dxMat = fvm.getDxMat();
  const SparseMatrix<double> &dyMat = fvm.getDyMat();

  REQUIRE(dxMat.rows() == n);
  REQUIRE(dyMat.rows() == n);

  VectorXd values(n);
  values << 1, 2, 3, 0;

  VectorXd dxs = dxMat * values;
  VectorXd dys = dyMat * values;

  REQUIRE(dxs(0) == 1);
  REQUIRE(dys(0) == 1);
}

TEST_CASE("FVM createGradientMatrix: Basic", "[Gradient]") {
  //          0 1 3
  //          - - -
  // 3 4 5    . . . | 3
  // 0 1 2    . . . | 1
  int n = 6;
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

  const SparseMatrix<double> &dxMat = fvm.getDxMat();
  const SparseMatrix<double> &dyMat = fvm.getDyMat();

  REQUIRE(dxMat.rows() == n);
  REQUIRE(dyMat.rows() == n);
}