#include <catch2/catch_test_macros.hpp>

#include "../src/FVM.h"
#include "../src/Grid.h"
#include "../src/Surface.h"
#include "../src/Debug.h"

#include "GridBuilder.h"

TEST_CASE("PDE - Poisson's Equation")
{
  int rows = 5;
  int cols = 5;
  double w = 1;
  double h = 1;
  vector<Face> faces = GridBuilder::buildRectangularGrid(rows, cols, w, h);

  Grid grid = Grid(faces);
  int n = grid.getN();
  int m = grid.getM();

  vector<std::pair<int, double>> fixed(m);
  for (int i = 0; i < grid.getM(); i++)
  {
    fixed[i] = {i, 0};
  }
  Surface surface = Surface(fixed);
  surface.init(faces, grid);

  FVM fvm(grid);

  VectorXd phi = VectorXd::Zero(n);
  VectorXd f = VectorXd::Ones(n);

  int iterations = 10;
  while (iterations--)
  {
    VectorXd gamma = VectorXd::Ones(n);

    auto [M, b] = fvm.laplacian(gamma, surface);

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
  vector<Face> faces = GridBuilder::buildRectangularGrid(rows, cols, w, h);

  Grid grid = Grid(faces);
  int n = grid.getN();
  int m = grid.getM();

  vector<std::pair<int, double>> fixed(m);
  for (int i = 0; i < m; i++)
  {
    fixed[i] = {i, 0};
  }
  Surface surface = Surface(fixed);
  surface.init(faces, grid);

  FVM fvm(grid);

  VectorXd phi = VectorXd::Zero(n);

  auto [M, b] = fvm.laplacian(VectorXd::Ones(n), surface);
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

TEST_CASE("PDE - Convergence: Transport Equation")
{
  int rows = 5;
  int cols = 5;
  double w = 1;
  double h = 1;
  vector<Face> faces = GridBuilder::buildRectangularGrid(rows, cols, w, h);

  Grid grid = Grid(faces);
  int n = grid.getN();
  int m = grid.getM();

  Surface::init(faces, grid);

  vector<std::pair<int, double>> pFixed;
  for (int i = 0; i < m; i++)
    pFixed.push_back({i, 0});
  Surface pSf = Surface(pFixed);

  vector<std::pair<int, double>> uFixed;
  for (int i = 0; i < m; i++)
    uFixed.push_back({i, 0});
  Surface uSf = Surface(uFixed);

  vector<std::pair<int, double>> vFixed;
  for (int i = 0; i < m; i++)
    vFixed.push_back({i, 0});
  Surface vSf = Surface(vFixed);

  FVM fvm(grid);

  VectorXd p = VectorXd::Ones(n);
  VectorXd u = VectorXd::Ones(n);
  VectorXd v = VectorXd::Ones(n);
  VectorXd s = VectorXd::Ones(n);

  auto [M, b] = fvm.laplacian(VectorXd::Ones(n), pSf);
  Eigen::BiCGSTAB<SpMat> solver;
  solver.compute(M);

  double alpha = 0.6;

  int iterations = 10;
  while (iterations--)
  {
    auto [pDx, pDy] = fvm.calcDf(p, pSf);
    VectorXd b = fvm.calcDiv(u, v, uSf, vSf).cwiseProduct(p) + pDx.cwiseProduct(u) + pDy.cwiseProduct(v) + s;

    VectorXd res = solver.solve(b);

    Debug::debugResiduals(res, p, "");

    p = (1 - alpha) * p + (alpha)*res;

    Debug::debug2d(p, rows, cols);
  }
}

TEST_CASE("PDE - Convergence: Velocity Transport")
{
  int rows = 5;
  int cols = 5;
  double w = 1;
  double h = 1;
  vector<Face> faces = GridBuilder::buildRectangularGrid(rows, cols, w, h);

  Grid grid = Grid(faces);
  int n = grid.getN();
  int m = grid.getM();

  Surface::init(faces, grid);

  vector<std::pair<int, double>> uFixed;
  for (int i = 0; i < m; i++)
    uFixed.push_back({i, 1});
  Surface uSf = Surface(uFixed);

  vector<std::pair<int, double>> vFixed;
  for (int i = 0; i < m; i++)
    vFixed.push_back({i, 0});
  Surface vSf = Surface(vFixed);

  FVM fvm(grid);

  VectorXd u = VectorXd::Zero(n);
  VectorXd v = VectorXd::Ones(n);

  auto [Mu, bu] = fvm.laplacian(VectorXd::Ones(n), uSf);
  Eigen::BiCGSTAB<SpMat> uSolver;
  uSolver.compute(Mu);

  auto [Mv, bv] = fvm.laplacian(VectorXd::Ones(n), vSf);
  Eigen::BiCGSTAB<SpMat> vSolver;
  vSolver.compute(Mv);

  double alpha = 0.6;

  int iterations = 10;
  while (iterations--)
  {
    auto [uDx, uDy] = fvm.calcDf(u, uSf);
    VectorXd fu = fvm.calcDiv(u, v, uSf, vSf).cwiseProduct(u) + uDx.cwiseProduct(u) + uDy.cwiseProduct(v);

    auto [vDx, vDy] = fvm.calcDf(v, vSf);
    VectorXd fv = fvm.calcDiv(u, v, uSf, vSf).cwiseProduct(v) + vDx.cwiseProduct(u) + vDy.cwiseProduct(v);

    VectorXd uRes = uSolver.solve(fu - bu);
    VectorXd vRes = vSolver.solve(fv - bv);

    Debug::debugResiduals(uRes, u, "u");
    Debug::debugResiduals(vRes, v, "v");

    u = (1 - alpha) * u + (alpha)*uRes;
    v = (1 - alpha) * v + (alpha)*vRes;

    Debug::debug2d(u, rows, cols);
    Debug::debug2d(v, rows, cols);
  }
}

TEST_CASE("PDE - Convergence: Velocity Transport (1)")
{
  int rows = 5;
  int cols = 5;
  double w = 1;
  double h = 1;
  vector<Face> faces = GridBuilder::buildRectangularGrid(rows, cols, w, h);

  Grid grid = Grid(faces);
  int n = grid.getN();
  int m = grid.getM();

  Surface::init(faces, grid);

  vector<std::pair<int, double>> uFixed;
  for (int i = 0; i < m; i++)
    uFixed.push_back({i, 1});
  Surface uSf = Surface(uFixed);

  vector<std::pair<int, double>> vFixed;
  for (int i = 0; i < m; i++)
    vFixed.push_back({i, -1});
  Surface vSf = Surface(vFixed);

  FVM fvm(grid);

  VectorXd u = VectorXd::Zero(n);
  VectorXd v = VectorXd::Ones(n);

  auto [Mu, bu] = fvm.laplacian(VectorXd::Ones(n), uSf);
  Eigen::BiCGSTAB<SpMat> uSolver;
  uSolver.compute(Mu);

  auto [Mv, bv] = fvm.laplacian(VectorXd::Ones(n), vSf);
  Eigen::BiCGSTAB<SpMat> vSolver;
  vSolver.compute(Mv);

  double alpha = 0.6;

  int iterations = 20;
  while (iterations--)
  {
    auto [uDx, uDy] = fvm.calcDf(u, uSf);
    VectorXd fu = fvm.calcDiv(u, v, uSf, vSf).cwiseProduct(u) + uDx.cwiseProduct(u) + uDy.cwiseProduct(v);

    auto [vDx, vDy] = fvm.calcDf(v, vSf);
    VectorXd fv = fvm.calcDiv(u, v, uSf, vSf).cwiseProduct(v) + vDx.cwiseProduct(u) + vDy.cwiseProduct(v);

    VectorXd uRes = uSolver.solve(fu - bu);
    VectorXd vRes = vSolver.solve(fv - bv);

    Debug::debugResiduals(uRes, u, "u");
    Debug::debugResiduals(vRes, v, "v");

    u = (1 - alpha) * u + (alpha)*uRes;
    v = (1 - alpha) * v + (alpha)*vRes;

    Debug::debug2d(u, rows, cols);
    Debug::debug2d(v, rows, cols);
  }
}

TEST_CASE("PDE - Convergence: Interdependence")
{
  int rows = 5;
  int cols = 5;
  double w = 1;
  double h = 1;
  vector<Face> faces = GridBuilder::buildRectangularGrid(rows, cols, w, h);

  Grid grid = Grid(faces);
  int n = grid.getN();
  int m = grid.getM();

  Surface::init(faces, grid);

  vector<std::pair<int, double>> pFixed;
  for (int i = 0; i < m; i++)
    pFixed.push_back({i, 1});
  Surface pSf = Surface(pFixed);

  vector<std::pair<int, double>> uFixed;
  for (int i = 0; i < m; i++)
  {
    if (faces[i].c.x() == 0)
      uFixed.push_back({i, 1});
    if (faces[i].c.x() == cols)
      uFixed.push_back({i, 2});
  }
  Surface uSf = Surface(uFixed);

  vector<std::pair<int, double>> vFixed;
  for (int i = 0; i < m; i++)
    vFixed.push_back({i, -1});
  Surface vSf = Surface(vFixed);

  FVM fvm(grid);

  VectorXd u = VectorXd::Zero(n);
  VectorXd v = VectorXd::Zero(n);
  VectorXd p = VectorXd::Zero(n);

  auto [Mu, bu] = fvm.laplacian(VectorXd::Ones(n), uSf);
  Eigen::BiCGSTAB<SpMat> uSolver;
  uSolver.compute(Mu);

  auto [Mv, bv] = fvm.laplacian(VectorXd::Ones(n), vSf);
  Eigen::BiCGSTAB<SpMat> vSolver;
  vSolver.compute(Mv);

  auto [Mp, bp] = fvm.laplacian(VectorXd::Ones(n), pSf);
  Eigen::BiCGSTAB<SpMat> pSolver;
  pSolver.compute(Mp);

  VectorXd uS(n);
  for (int i = 0; i < rows; i++)
  {
    uS(i * cols + 3) = 1;
  }

  double alpha = 0.6;

  int iterations = 20;
  while (iterations--)
  {
    auto [uDx, uDy] = fvm.calcDf(u, uSf);
    VectorXd fu = fvm.calcDiv(u, v, uSf, vSf).cwiseProduct(u) + uDx.cwiseProduct(u) + uDy.cwiseProduct(v) + uS;

    auto [vDx, vDy] = fvm.calcDf(v, vSf);
    VectorXd fv = fvm.calcDiv(u, v, uSf, vSf).cwiseProduct(v) + vDx.cwiseProduct(u) + vDy.cwiseProduct(v);

    VectorXd uRes = uSolver.solve(fu - bu);
    VectorXd vRes = vSolver.solve(fv - bv);

    Debug::debugResiduals(uRes, u, "u");
    Debug::debugResiduals(vRes, v, "v");

    u = (1 - alpha) * u + (alpha)*uRes;
    v = (1 - alpha) * v + (alpha)*vRes;

    auto [pDx, pDy] = fvm.calcDf(p, pSf);
    VectorXd fp = fvm.calcDiv(u, v, uSf, vSf).cwiseProduct(p) + pDx.cwiseProduct(u) + pDy.cwiseProduct(v);
    VectorXd pRes = pSolver.solve(fp - bp);

    Debug::debugResiduals(pRes, p, "p");

    p = (1 - alpha) * p + (alpha)*pRes;

    Debug::debug2d(u, rows, cols, "u");
    Debug::debug2d(v, rows, cols, "v");
    Debug::debug2d(p, rows, cols, "p");
  }
}