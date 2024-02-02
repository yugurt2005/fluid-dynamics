#include <catch2/catch_test_macros.hpp>

#include "../src/FVM.h"
#include "../src/Grid.h"
#include "../src/Surface.h"
#include "../src/Debug.h"

#include "GridBuilder.h"

TEST_CASE("FVM - laplacian: f(x) = x*x - 2", "[laplacian]")
{
  int n = 1;
  int m = 3;
  double w = 1;
  double h = 1;
  vector<Face> faces = GridBuilder::buildRectangularGrid(n, m, w, h);

  Grid grid = Grid(faces);

  vector<std::pair<int, double>> fixed(grid.getM());
  for (int i = 0; i < grid.getM(); i++)
  {
    fixed[i] = {i, 0};
  }
  Surface surface = Surface(fixed);
  surface.init(faces, grid);

  VectorXd gamma = VectorXd::Ones(grid.getM());

  VectorXd phi(n * m);
  phi << -2, 0, 6;

  auto [res, _] = FVM(grid).laplacian(gamma, surface);
  VectorXd ans = res * phi;

  REQUIRE(res.rows() == n * m);
  REQUIRE(res.cols() == n * m);

  CHECK(ans(1) == 4);
}

TEST_CASE("FVM - laplacian: f(x, y) = x + y", "[laplacian]")
{
  int n = 4;
  int m = 4;
  double w = 1;
  double h = 1;
  vector<Face> faces = GridBuilder::buildRectangularGrid(n, m, w, h);

  Grid grid = Grid(faces);

  vector<std::pair<int, double>> fixed(grid.getM());
  for (int i = 0; i < grid.getM(); i++)
  {
    fixed[i] = {i, 0};
  }
  Surface surface = Surface(fixed);
  surface.init(faces, grid);

  VectorXd gamma = VectorXd::Ones(grid.getM());

  VectorXd phi(n * m);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      phi(i * m + j) = i + j + 2;
    }
  }

  auto [res, _] = FVM(grid).laplacian(gamma, surface);
  VectorXd ans = res * phi;

  for (int i = 1; i < n - 1; i++)
  {
    for (int j = 1; j < m - 1; j++)
    {
      CHECK(ans(i * m + j) == 0);
    }
  }
}

TEST_CASE("FVM - laplacian: f(x, y) = x*x * y*y", "[laplacian]")
{
  auto f = [](double i, double j)
  {
    return i * i * j * j;
  };

  SECTION("Regular")
  {
    int n = 3;
    int m = 4;
    double w = 1;
    double h = 1;
    vector<Face> faces = GridBuilder::buildRectangularGrid(n, m, w, h);

    Grid grid = Grid(faces);

    vector<std::pair<int, double>> fixed(grid.getM());
    for (int i = 0; i < grid.getM(); i++)
    {
      fixed[i] = {i, 0};
    }
    Surface surface = Surface(fixed);
    surface.init(faces, grid);

    VectorXd gamma = VectorXd::Ones(grid.getM());

    VectorXd phi(n * m);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < m; j++)
        phi(i * m + j) = f(i, j);

    auto [res, _] = FVM(grid).laplacian(gamma, surface);
    VectorXd ans = res * phi;

    for (int i = 1; i < n - 1; i++)
    {
      for (int j = 1; j < m - 1; j++)
      {
        double x = i;
        double y = j;
        CHECK(ans(i * m + j) == 2 * y * y + 2 * x * x);
      }
    }
  }

  SECTION("Small Scale")
  {
    double scale = 0.1;

    int n = 3;
    int m = 4;
    double w = scale;
    double h = scale;
    vector<Face> faces = GridBuilder::buildRectangularGrid(n, m, w, h);

    Grid grid = Grid(faces);

    vector<std::pair<int, double>> fixed(grid.getM());
    for (int i = 0; i < grid.getM(); i++)
    {
      fixed[i] = {i, 0};
    }
    Surface surface = Surface(fixed);
    surface.init(faces, grid);

    VectorXd gamma = VectorXd::Ones(grid.getM());

    VectorXd phi(n * m);
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < m; j++)
      {
        phi(i * m + j) = f(i * scale, j * scale);
      }
    }

    auto [res, _] = FVM(grid).laplacian(gamma, surface);
    VectorXd ans = res * phi;

    for (int i = 1; i < n - 1; i++)
    {
      for (int j = 1; j < m - 1; j++)
      {
        double x = i * scale;
        double y = j * scale;
        CHECK(std::abs(ans(i * m + j) - (2 * y * y + 2 * x * x) * grid.getVolume(i * m + j)) < 1e-9);
      }
    }
  }
}

TEST_CASE("FVM - laplacian: Mass Flux", "[laplacian]") {
  int n = 3;
  int m = 4;
  double w = 1;
  double h = 1;
  vector<Face> faces = GridBuilder::buildRectangularGrid(n, m, w, h);

  Grid grid = Grid(faces);

  vector<std::pair<int, double>> fixed;
  for (int i = 0; i < faces.size(); i++) {
    double x = faces[i].c.x();
    double y = faces[i].c.y();
    if (y == 0 || y == n || x == 0 || x == m) {
      fixed.push_back({i, x + y});
    }
  }

  Surface surface = Surface(fixed);
  surface.init(faces, grid);

  VectorXd gamma = VectorXd::Ones(grid.getM());

  VectorXd phi(n * m);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      phi(i * m + j) = i + j + 1;

  auto [M, b] = FVM(grid).laplacian(gamma, surface);
  VectorXd ans = M * phi + b;

  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      CHECK(ans(i * m + j) == 0);
    }
  }
}