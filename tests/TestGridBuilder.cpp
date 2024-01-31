#include <catch2/catch_test_macros.hpp>

#include "GridBuilder.h"

TEST_CASE("GridBuilder - buildRectangularGrid: Happy Test") {
  int n = 2;
  int m = 2;
  double w = 1;
  double h = 1;
  Grid grid = GridBuilder::buildRectangularGrid(n, m, w, h);

  REQUIRE(grid.getN() == n * m);
  REQUIRE(grid.getM() == n * (m + 1) + m * (n + 1));
}

TEST_CASE("GridBuilder - buildRectangularGrid: Cell Centers") {
  SECTION("1 x 1") {
    int n = 2;
    int m = 2;
    double w = 1;
    double h = 1;
    Grid grid = GridBuilder::buildRectangularGrid(n, m, w, h);

    CHECK(grid.getCenter(0) == Vector2d(0.5, 0.5));
  }

  SECTION("2.5 x 2.5") {
    int n = 2;
    int m = 2;
    double w = 2.5;
    double h = 2.5;
    Grid grid = GridBuilder::buildRectangularGrid(n, m, w, h);

    CHECK(grid.getCenter(0) == Vector2d(1.25, 1.25));
  }

  SECTION("1.5 x 2.5") {
    int n = 2;
    int m = 2;
    double w = 1.5;
    double h = 2.5;
    Grid grid = GridBuilder::buildRectangularGrid(n, m, w, h);

    CHECK(grid.getCenter(0) == Vector2d(0.75, 1.25));
  }
}

TEST_CASE("GridBuilder - buildRectangularGrid: Face Normals") {
  SECTION("1.5 x 2.5") {
    int n = 2;
    int m = 2;
    double w = 1.5;
    double h = 2.5;
    Grid grid = GridBuilder::buildRectangularGrid(n, m, w, h);

    CHECK(grid.getFaces()[0].n == Vector2d(0, -1));
  }
}

TEST_CASE("GridBuilder - buildRectangularGrid: Edge Normals") {
  SECTION("1.5 x 2.5") {
    int n = 2;
    int m = 2;
    double w = 2.5;
    double h = 2.5;
    Grid grid = GridBuilder::buildRectangularGrid(n, m, w, h);

    CHECK(grid.getAdj(0)[1].n == Vector2d(0, 1));
    CHECK(grid.getAdj(2)[0].n == Vector2d(0, -1));
  }
}

TEST_CASE("GridBuilder - buildTriangularGrid: Happy Test") {
  int n = 2;
  int m = 2;
  double w = 1;
  double h = 1;
  Grid grid = GridBuilder::buildTriangularGrid(n, m, w, h);

  REQUIRE(grid.getN() == n * m * 2);
  REQUIRE(grid.getM() == n * (m + 1) + m * (n + 1) + n * m);
}

TEST_CASE("GridBuilder - buildTriangularGrid: Cell Centers") {
  SECTION("1 x 1") {
    int n = 2;
    int m = 2;
    double w = 1;
    double h = 1;
    Grid grid = GridBuilder::buildTriangularGrid(n, m, w, h);

    CHECK(grid.getCenter(0) == Vector2d(1.0/3, 1.0/3));
  }

  SECTION("2.5 x 2.5") {
    int n = 2;
    int m = 2;
    double w = 2.5;
    double h = 2.5;
    Grid grid = GridBuilder::buildTriangularGrid(n, m, w, h);

    CHECK(grid.getCenter(0) == Vector2d(2.5/3, 2.5/3));
  }

  SECTION("1.5 x 2.5") {
    int n = 2;
    int m = 2;
    double w = 1.5;
    double h = 2.5;
    Grid grid = GridBuilder::buildTriangularGrid(n, m, w, h);

    CHECK(grid.getCenter(0) == Vector2d(1.5/3, 2.5/3));
  }
}

TEST_CASE("GridBuilder - buildTriangularGrid: Face Normals") {
  SECTION("1 x 2") {
    int n = 2;
    int m = 2;
    double w = 1;
    double h = 1;
    Grid grid = GridBuilder::buildTriangularGrid(n, m, w, h);

    CHECK(grid.getFaces()[0].n == Vector2d(1, 1).normalized());
  }

  SECTION("1.5 x 2.5") {
    int n = 2;
    int m = 2;
    double w = 1.5;
    double h = 2.5;
    Grid grid = GridBuilder::buildTriangularGrid(n, m, w, h);

    CHECK(grid.getFaces()[0].n == Vector2d(2.5, 1.5).normalized());
  }
}

TEST_CASE("GridBuilder - buildTriangularGrid: Edge Normals") {
  SECTION("1.5 x 2.5") {
    int n = 2;
    int m = 2;
    double w = 2.5;
    double h = 2.5;
    Grid grid = GridBuilder::buildTriangularGrid(n, m, w, h);

    CHECK(grid.getAdj(0)[0].n == Vector2d(1, 1).normalized());
    CHECK(grid.getAdj(1)[0].n == Vector2d(-1, -1).normalized());
  }
}

// TODO: distance tests