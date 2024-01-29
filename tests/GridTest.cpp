#include <catch2/catch_test_macros.hpp>

#include <iostream>

#include "GridHelper.h"

TEST_CASE("Grid Helper")
{
  SECTION("Build Rectangular Grid: 1 x 1")
  {
    Grid grid = buildRectangularGrid(1, 1, 1);

    CHECK(grid.getN() == 1);
    CHECK(grid.getZ() == 4);
  }

  SECTION("Build Rectangular Grid: 2 x 1")
  {
    Grid grid = buildRectangularGrid(2, 1, 1);

    CHECK(grid.getN() == 2);
    CHECK(grid.getZ() == 7);
  }

  SECTION("Build Rectangular Grid: 2 x 2")
  {
    Grid grid = buildRectangularGrid(2, 2, 1);

    CHECK(grid.getN() == 4);
    CHECK(grid.getZ() == 12);
  }

  SECTION("Build Rectangular Grid: 2 x 3")
  {
    int n = 2;
    int m = 3;
    Grid grid = buildRectangularGrid(n, m, 1);

    CHECK(grid.getN() == n * m);
    CHECK(grid.getZ() == n * (m + 1) + m * (n + 1));
  }

  SECTION("Build Rectangular Grid: 3 x 3")
  {
    int n = 3;
    int m = 3;
    Grid grid = buildRectangularGrid(n, m, 1);

    CHECK(grid.getN() == n * m);
    CHECK(grid.getZ() == n * (m + 1) + m * (n + 1));
  }
}

TEST_CASE("Grid: Centers")
{
  SECTION("Rectangular Grid: 2 x 1")
  {
    Grid grid = buildRectangularGrid(2, 1, 1);

    CHECK(grid.getCenter(0) == Vector2d(0.5, 0.5));
    CHECK(grid.getCenter(1) == Vector2d(0.5, 1.5));
  }

  SECTION("Rectangular Grid: 1 x 2")
  {
    Grid grid = buildRectangularGrid(1, 2, 1);

    CHECK(grid.getCenter(0) == Vector2d(0.5, 0.5));
    CHECK(grid.getCenter(1) == Vector2d(1.5, 0.5));
  }

  SECTION("Rectangular Grid: 3 x 3")
  {
    Grid grid = buildRectangularGrid(3, 3, 1);

    CHECK(grid.getCenter(2) == Vector2d(2.5, 0.5));
    CHECK(grid.getCenter(3) == Vector2d(0.5, 1.5));
    CHECK(grid.getCenter(8) == Vector2d(2.5, 2.5));
  }
}

TEST_CASE("Grid: Neighbors")
{
  SECTION("Unique")
  {
    Grid grid = buildRectangularGrid(3, 3, 1);

    int n = 9;
    for (int c = 0; c < n; c++)
    {
      vector<int> neighbors = grid.getNeighbors(c);
      for (int i = 0; i < neighbors.size(); i++)
      {
        for (int j = 0; j < neighbors.size(); j++)
        {
          if (i != j)
            CHECK(neighbors[i] != neighbors[j]);
        }
      }
    }
  }
}
