#include <catch2/catch_test_macros.hpp>

#include "../src/FVM.h"
#include "../src/Debug.h"

#include "GridBuilder.h"

TEST_CASE("PDE - Poisson's Equation")
{
  int rows = 5;
  int cols = 5;
  double w = 1;
  double h = 1;
  Grid grid = GridBuilder::buildRectangularGrid(rows, cols, w, h);

  FVM fvm(grid);

  int n = grid.getN();
  int m = grid.getM();

  VectorXd phi = VectorXd::Zero(n);
  VectorXd f = VectorXd::Ones(n);

  int iterations = 10;
  while (iterations--)
  {
    VectorXd gamma = VectorXd::Ones(n);

    SpMat M = fvm.laplacian(gamma);

    Eigen::BiCGSTAB<SpMat> solver;
    solver.compute(M);

    phi = solver.solve(f);

    Debug::debug2d(phi, rows, cols);
  }
}

TEST_CASE("PDE - Convergence: Iterative Solution for One Variable")
{
  int rows = 5;
  int cols = 5;
  double w = 1;
  double h = 1;
  Grid grid = GridBuilder::buildRectangularGrid(rows, cols, w, h);

  FVM fvm(grid);

  int n = grid.getN();
  int m = grid.getM();

  VectorXd phi = VectorXd::Zero(n);

  SpMat M = fvm.laplacian(VectorXd::Ones(n));
  Eigen::BiCGSTAB<SpMat> solver;
  solver.compute(M);

  double alpha = 0.6;

  int iterations = 10;
  while (iterations--)
  {
    VectorXd f = phi + VectorXd::Ones(n);
    VectorXd res = solver.solve(f);

    Debug::debugResiduals(res, phi, "");

    phi = (1 - alpha) * phi + (alpha)*res;

    Debug::debug2d(phi, rows, cols);
  }

  VectorXd lhs = M * phi;
  VectorXd rhs = phi + VectorXd::Ones(n);

  CHECK(Debug::debugResiduals(lhs, rhs, "Total Error") < 1e-3);
}