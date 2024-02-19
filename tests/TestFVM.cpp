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

TEST_CASE("FVM - laplacian: Boundaries", "[laplacian]")
{
  int n = 3;
  int m = 4;
  double w = 1;
  double h = 1;
  vector<Face> faces = GridBuilder::buildRectangularGrid(n, m, w, h);

  Grid grid = Grid(faces);

  vector<std::pair<int, double>> fixed;
  for (int i = 0; i < faces.size(); i++)
  {
    double x = faces[i].c.x();
    double y = faces[i].c.y();
    if (y == 0 || y == n || x == 0 || x == m)
    {
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

TEST_CASE("FVM - gradient: x + y", "[gradient]")
{
  int n = 3;
  int m = 4;
  double w = 1;
  double h = 1;
  vector<Face> faces = GridBuilder::buildRectangularGrid(n, m, w, h);

  Grid grid = Grid(faces);

  auto f = [](double i, double j)
  {
    return i + j;
  };

  vector<std::pair<int, double>> fixed;
  for (int i = 0; i < faces.size(); i++)
  {
    double x = faces[i].c.x();
    double y = faces[i].c.y();
    if (y == 0 || y == n || x == 0 || x == m)
    {
      fixed.push_back({i, f(x, y)});
    }
  }
  Surface surface = Surface(fixed);
  surface.init(faces, grid);

  VectorXd phi(n * m);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      phi(i * m + j) = f(i + 0.5, j + 0.5);

  FVM fvm = FVM(grid);

  auto [dx, dy] = fvm.calcDf(phi, surface);

  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      CHECK(dx(i * m + j) == 1);
      CHECK(dy(i * m + j) == 1);
    }
  }
}

TEST_CASE("FVM - gradient: x*x + y*y", "[gradient]")
{
  int n = 3;
  int m = 4;
  double w = 1;
  double h = 1;
  vector<Face> faces = GridBuilder::buildRectangularGrid(n, m, w, h);

  Grid grid = Grid(faces);

  auto f = [](Vector position)
  {
    double i = position.x();
    double j = position.y();
    return i * i + j * j;
  };

  vector<std::pair<int, double>> fixed;
  for (int i = 0; i < faces.size(); i++)
  {
    if (faces[i].isWall())
    {
      fixed.push_back({i, f(faces[i].c)});
    }
  }
  Surface surface = Surface(fixed);
  surface.init(faces, grid);

  VectorXd phi(n * m);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      phi(i * m + j) = f(grid.getCenter(i * m + j));

  FVM fvm = FVM(grid);

  auto [dx, dy] = fvm.calcDf(phi, surface);

  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      CHECK(std::abs(dx(i * m + j) - 2 * (j + 0.5)) <= 0.25);
      CHECK(std::abs(dy(i * m + j) - 2 * (i + 0.5)) <= 0.25);
    }
  }
}

TEST_CASE("FVM - gradient: x*y", "[gradient]")
{
  int n = 3;
  int m = 4;
  double w = 1;
  double h = 1;
  vector<Face> faces = GridBuilder::buildRectangularGrid(n, m, w, h);

  Grid grid = Grid(faces);

  auto f = [](Vector position)
  {
    double i = position.x();
    double j = position.y();
    return i * j;
  };

  vector<std::pair<int, double>> fixed;
  for (int i = 0; i < faces.size(); i++)
  {
    if (faces[i].isWall())
    {
      fixed.push_back({i, f(faces[i].c)});
    }
  }
  Surface surface = Surface(fixed);
  surface.init(faces, grid);

  VectorXd phi(n * m);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      phi(i * m + j) = f(grid.getCenter(i * m + j));

  FVM fvm = FVM(grid);

  auto [dx, dy] = fvm.calcDf(phi, surface);

  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      CHECK(dx(i * m + j) == i + 0.5);
      CHECK(dy(i * m + j) == j + 0.5);
    }
  }
}

TEST_CASE("FVM - div #1", "[div]")
{
  int n = 3;
  int m = 4;
  double w = 1;
  double h = 1;
  vector<Face> faces = GridBuilder::buildRectangularGrid(n, m, w, h);

  Grid grid = Grid(faces);

  auto p = [](Vector position)
  {
    double i = position.x();
    double j = position.y();
    return i * j;
  };

  auto q = [](Vector position)
  {
    double i = position.x();
    double j = position.y();
    return i + j;
  };

  auto f = [](Vector position)
  {
    double i = position.x();
    double j = position.y();
    return j + 1;
  };

  vector<std::pair<int, double>> uFix;
  for (int i = 0; i < faces.size(); i++)
  {
    if (faces[i].isWall())
    {
      uFix.push_back({i, p(faces[i].c)});
    }
  }

  vector<std::pair<int, double>> vFix;
  for (int i = 0; i < faces.size(); i++)
  {
    if (faces[i].isWall())
    {
      vFix.push_back({i, q(faces[i].c)});
    }
  }
  Surface::init(faces, grid);
  Surface uSf(uFix);
  Surface vSf(vFix);

  VectorXd u(n * m);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      u(i * m + j) = p(grid.getCenter(i * m + j));
  VectorXd v(n * m);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      v(i * m + j) = q(grid.getCenter(i * m + j));

  Debug::debug2d(u, n, m);
  Debug::debug2d(v, n, m);

  FVM fvm = FVM(grid);

  VectorXd res = fvm.calcDiv(u, v, uSf, vSf);

  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      CHECK(res(i * m + j) == f(grid.getCenter(i * m + j)));
    }
  }
}

TEST_CASE("FVM - div: #2", "[div]")
{
  int n = 3;
  int m = 4;
  double w = 1;
  double h = 1;
  vector<Face> faces = GridBuilder::buildRectangularGrid(n, m, w, h);

  Grid grid = Grid(faces);

  auto p = [](Vector position)
  {
    double x = position.x();
    double y = position.y();
    return x * x + y * y;
  };

  auto q = [](Vector position)
  {
    double x = position.x();
    double y = position.y();
    return 0.5 * x + 2 * y / (x + 10);
  };

  auto f = [](Vector position)
  {
    double x = position.x();
    double y = position.y();
    return 2 * x + 2 / (x + 10);
  };

  vector<std::pair<int, double>> uFix;
  for (int i = 0; i < faces.size(); i++)
  {
    if (faces[i].isWall())
    {
      uFix.push_back({i, p(faces[i].c)});
    }
  }

  vector<std::pair<int, double>> vFix;
  for (int i = 0; i < faces.size(); i++)
  {
    if (faces[i].isWall())
    {
      vFix.push_back({i, q(faces[i].c)});
    }
  }

  Surface::init(faces, grid);
  Surface uSf(uFix);
  Surface vSf(vFix);

  VectorXd u(n * m);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      u(i * m + j) = p(grid.getCenter(i * m + j));

  VectorXd v(n * m);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      v(i * m + j) = q(grid.getCenter(i * m + j));

  Debug::debug2d(u, n, m);
  Debug::debug2d(v, n, m);

  FVM fvm = FVM(grid);

  VectorXd res = fvm.calcDiv(u, v, uSf, vSf);

  Debug::debug2d(res, n, m);

  assert(res.size() == n * m);

  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      double difference = res(i * m + j) - f(grid.getCenter(i * m + j));

      CHECK(std::abs(difference) <= 0.3);
    }
  }
}

TEST_CASE("FVM - buildGradients: x*y") {
  int n = 3;
  int m = 4;
  double w = 1;
  double h = 1;
  vector<Face> faces = GridBuilder::buildRectangularGrid(n, m, w, h);

  Grid grid = Grid(faces);

  auto f = [](Vector position)
  {
    double x = position.x();
    double y = position.y();
    return x * y;
  };
  auto fx = [](Vector position) {
    return position.y();
  };
  auto fy = [](Vector position) {
    return position.x();
  };

  vector<std::pair<int, double>> fixed;
  for (int i = 0; i < faces.size(); i++)
  {
    if (faces[i].isWall())
    {
      fixed.push_back({i, f(faces[i].c)});
    }
  }
  Surface surface = Surface(fixed);
  surface.init(faces, grid);

  VectorXd phi(n * m);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      phi(i * m + j) = f(grid.getCenter(i * m + j));

  FVM fvm = FVM(grid);

  auto [Gx, bx, Gy, by] = fvm.buildGradients(surface);

  VectorXd dx = Gx * phi + bx;
  VectorXd dy = Gy * phi + by;

  // Debug::debug2d(phi, n, m, "phi");
  // Debug::debug2d(dx, n, m, "dx");
  // Debug::debug2d(dy, n, m, "dy");

  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      CHECK(dx(i * m + j) == fx(grid.getCenter(i * m + j)));
      CHECK(dy(i * m + j) == fy(grid.getCenter(i * m + j)));
    }
  } 
}

TEST_CASE("FVM - buildGradients: x + y") {
  int n = 3;
  int m = 4;
  double w = 1;
  double h = 1;
  vector<Face> faces = GridBuilder::buildRectangularGrid(n, m, w, h);

  Grid grid = Grid(faces);

  auto f = [](Vector position)
  {
    double x = position.x();
    double y = position.y();
    return x + y;
  };
  auto fx = [](Vector position) {
    return 1;
  };
  auto fy = [](Vector position) {
    return 1;
  };

  vector<std::pair<int, double>> fixed;
  for (int i = 0; i < faces.size(); i++)
  {
    if (faces[i].isWall())
    {
      fixed.push_back({i, f(faces[i].c)});
    }
  }
  Surface surface = Surface(fixed);
  surface.init(faces, grid);

  VectorXd phi(n * m);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      phi(i * m + j) = f(grid.getCenter(i * m + j));

  FVM fvm = FVM(grid);

  auto [Gx, bx, Gy, by] = fvm.buildGradients(surface);

  VectorXd dx = Gx * phi + bx;
  VectorXd dy = Gy * phi + by;

  // Debug::debug2d(phi, n, m, "phi");
  // Debug::debug2d(dx, n, m, "dx");
  // Debug::debug2d(dy, n, m, "dy");

  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      CHECK(dx(i * m + j) == fx(grid.getCenter(i * m + j)));
      CHECK(dy(i * m + j) == fy(grid.getCenter(i * m + j)));
    }
  } 
}

TEST_CASE("FVM - buildGradients: 3*x*x + 2*y*y") {
  // TODO: figure out boundaries

  int n = 4;
  int m = 4;
  double w = 1;
  double h = 1;
  vector<Face> faces = GridBuilder::buildRectangularGrid(n, m, w, h);

  Grid grid = Grid(faces);

  auto f = [](Vector position)
  {
    double x = position.x();
    double y = position.y();
    return 3 * x * x + 2 * y * y;
  };
  auto fx = [](Vector position) {
    double x = position.x();
    double y = position.y();
    return 6 * x;
  };
  auto fy = [](Vector position) {
    double x = position.x();
    double y = position.y();
    return 4 * y;
  };

  vector<std::pair<int, double>> fixed;
  for (int i = 0; i < faces.size(); i++)
  {
    if (faces[i].isWall())
    {
      fixed.push_back({i, f(faces[i].c)});
    }
  }
  Surface surface = Surface(fixed);
  surface.init(faces, grid);

  VectorXd phi(n * m);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      phi(i * m + j) = f(grid.getCenter(i * m + j));

  FVM fvm = FVM(grid);

  auto [Gx, bx, Gy, by] = fvm.buildGradients(surface);

  Debug::debug2d(bx, 1, n * m, "bx");
  Debug::debug2d(by, 1, n * m, "by");

  VectorXd dx = Gx * phi + bx;
  VectorXd dy = Gy * phi + by;

  VectorXd expectedDx = VectorXd::Zero(n * m);
  VectorXd expectedDy = VectorXd::Zero(n * m);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      expectedDx(i * m + j) = fx(grid.getCenter(i * m + j));
      expectedDy(i * m + j) = fy(grid.getCenter(i * m + j));
    }
  }

  Debug::debug2d(phi, n, m, "phi");
  Debug::debug2d(dx, n, m, "dx");
  Debug::debug2d(dy, n, m, "dy");
  Debug::debug2d(expectedDx, n, m, "expectedDx");
  Debug::debug2d(expectedDy, n, m, "expectedDy");

  for (int i = 1; i < n - 1; i++)
  {
    for (int j = 1; j < m - 1; j++)
    {
      CHECK(dx(i * m + j) == fx(grid.getCenter(i * m + j)));
      CHECK(dy(i * m + j) == fy(grid.getCenter(i * m + j)));
    }
  } 
}

TEST_CASE("FVM - calcMassFlux: u, v = x, y") {
  // TODO: figure out boundaries

  int n = 2;
  int m = 2;
  double w = 1;
  double h = 1;
  vector<Face> faces = GridBuilder::buildRectangularGrid(n, m, w, h);

  Grid grid = Grid(faces);

  auto uFun = [](Vector position)
  {
    double x = position.x();
    return x;
  };

  auto vFun = [](Vector position)
  {
    double y = position.y();
    return y;
  };

  Surface::init(faces, grid);

  vector<std::pair<int, double>> uFix;
  for (int i = 0; i < faces.size(); i++)
  {
    if (faces[i].isWall())
    {
      uFix.push_back({i, uFun(faces[i].c)});
    }
  }
  Surface uSf = Surface(uFix);

  vector<std::pair<int, double>> vFix;
  for (int i = 0; i < faces.size(); i++)
  {
    if (faces[i].isWall())
    {
      vFix.push_back({i, vFun(faces[i].c)});
    }
  }
  Surface vSf = Surface(vFix);

  VectorXd u(n * m);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      u(i * m + j) = uFun(grid.getCenter(i * m + j));

  VectorXd v(n * m);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      v(i * m + j) = vFun(grid.getCenter(i * m + j));

  FVM fvm = FVM(grid);

  VectorXd flux = fvm.calcMassFlux(u, v, uSf, vSf);

  for (int i = 0; i < grid.getM(); i++) {
    Face face = uSf.getFaces()[i];
    double expected = face.c.dot(face.n) * face.a;

    CHECK(flux(i) == expected);
    
    // cout << "flux(" << i << ") = " << flux(i) << " expected = " << expected << endl;
  }
}