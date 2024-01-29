#include <catch2/catch_test_macros.hpp>

#include "GridHelper.h"

#include "../src/FVM.h"

TEST_CASE("FVM Gradient: Happy Test")
{
  SECTION("Happy Test")
  {
    Grid grid = buildRectangularGrid(2, 2, 1);
    FVM fvm(grid);

    const SpMat &Gx = fvm.getGx();
    const SpMat &Gy = fvm.getGy();

    /*
    3 4
    1 2
    */
    VectorXd phi(4);
    phi << 1, 2, 3, 4;

    VectorXd dxs = Gx * phi;
    VectorXd dys = Gy * phi;

    CHECK(dxs(0) == 1);
    CHECK(dys(0) == 2);

    CHECK(dxs(2) == 1);
    CHECK(dys(2) == 2);
  }

  SECTION("Non-Uniform Phi")
  {
    Grid grid = buildRectangularGrid(3, 3, 1);
    FVM fvm(grid);

    const SpMat &Gx = fvm.getGx();
    const SpMat &Gy = fvm.getGy();

    /*
    5 4 6
    2 1 8
    5 3 9
    */
    VectorXd phi(9);
    phi << 5, 3, 9, 2, 1, 8, 5, 4, 6;

    VectorXd dxs = Gx * phi;
    VectorXd dys = Gy * phi;

    CHECK(dxs(0) == -2);
    CHECK(dys(0) == -3);

    CHECK(dxs(1) == 2);
    CHECK(dys(1) == -2);

    CHECK(dxs(4) == 3);
    CHECK(dys(4) == 0.5);
  }
}

TEST_CASE("FVM calcMassFlux: Happy Test")
{
  SECTION("Happy Test")
  {
    Grid grid = buildRectangularGrid(2, 2, 1);
    FVM fvm(grid);

    int n = grid.getN();
    int z = grid.getZ();

    VectorXd u(n);
    VectorXd v(n);
    u << 1, 1, 1, 1;
    v << 0, 0, 0, 0;

    VectorXd flux = fvm.calcMassFlux(u, v);

    CHECK(flux(0) == 1);
    CHECK(flux(1) == 0);
    CHECK(flux(2) == 0);
  }

  SECTION("Cancelling Fluxes")
  {
    Grid grid = buildRectangularGrid(2, 2, 1);
    FVM fvm(grid);

    int n = grid.getN();
    int z = grid.getZ();

    VectorXd u(n);
    VectorXd v(n);
    u << 1, -1, 1, -1;
    v << 1, 1, -1, -1;

    VectorXd flux = fvm.calcMassFlux(u, v);

    for (int i = 0; i < z; i++)
      CHECK(flux(i) == 0);
  }

  SECTION("div(Φ)")
  {
    Grid grid = buildRectangularGrid(2, 2, 1);
    FVM fvm(grid);

    int n = grid.getN();
    int z = grid.getZ();

    VectorXd u(n);
    VectorXd v(n);
    u << 2, 1, 2, 1;
    v << 0, 0, 0, 0;

    VectorXd flux = fvm.calcMassFlux(u, v);
    VectorXd div = fvm.getAdj() * flux;

    CHECK(flux(0) == 1.5);
    CHECK(flux(1) == 0);
    CHECK(flux(5) == 0);

    CHECK(div(0) == +1.5);
    CHECK(div(1) == -1.5);
  }

  SECTION("div(ΦΦ)") {
    // TODO
  }
}
