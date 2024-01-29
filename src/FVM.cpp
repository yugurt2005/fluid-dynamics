#include "FVM.h"

FVM::FVM(IGrid &grid) : grid(grid)
{
  n = grid.getN();
  z = grid.getZ();

  buildAdj();
  buildGradients();
}

SpMat FVM::convertDiagonal(const VectorXd &input)
{
  return input.asDiagonal().diagonal().sparseView();
}

void FVM::buildAdj()
{
  Adj = SpMat(n, z);

  vector<Face> faces = grid.getFaces();
  for (int i = 0; i < z; i++)
  {
    Face f = faces[i];

    if (f.l != -1)
      Adj.insert(f.l, i) = +1;
    if (f.r != -1)
      Adj.insert(f.r, i) = -1;
  }
}

void FVM::buildGradients()
{
  Gx = SpMat(n, n);
  Gy = SpMat(n, n);

  for (int c = 0; c < n; c++)
  {
    const vector<int> &neighbors = getNeighbors(c);

    vector<Vector2d> distances;
    for (int a : neighbors)
    {
      distances.push_back(grid.getCenter(a) - grid.getCenter(c));
    }

    for (int i = 0; i < neighbors.size(); i++)
    {
      MatrixXd gradient = calcLeastSquaresGradient(distances);

      Gx.insert(c, neighbors[i]) = gradient(0, i);
      Gy.insert(c, neighbors[i]) = gradient(1, i);

      Gx.coeffRef(c, c) -= gradient(0, i);
      Gy.coeffRef(c, c) -= gradient(1, i);
    }
  }
}

MatrixXd FVM::calcLeastSquaresGradient(const vector<Vector2d> &distances)
{
  int n = distances.size();

  MatrixXd d(n, 2);
  for (int i = 0; i < n; i++)
  {
    d(i, 0) = distances[i].x();
    d(i, 1) = distances[i].y();
  }

  MatrixXd w = MatrixXd::Zero(n, n);
  for (int i = 0; i < n; i++)
  {
    w(i, i) = 1 / distances[i].norm();
  }

  MatrixXd dT = d.transpose();
  MatrixXd wT = w.transpose();

  return (dT * wT * w * d).inverse() * dT * wT * w;
}

VectorXd FVM::calcMassFlux(const VectorXd &u, const VectorXd &v)
{
  assert(u.size() == n);
  assert(v.size() == n);

  vector<Face> faces = grid.getFaces();

  assert(faces.size() == z);

  VectorXd du = Gx * u;
  VectorXd dv = Gy * v;

  VectorXd uF(z);
  VectorXd vF(z);
  for (int i = 0; i < z; i++)
  {
    Face f = faces[i];

    // Upwind Differencing
    auto interpolate = [&](int index)
    {
      Vector2d center = grid.getCenter(index);
      double nu = u(index) + du(index) * (f.center.x() - center.x());
      double nv = v(index) + dv(index) * (f.center.y() - center.y());
      return Vector2d(nu, nv);
    };

    if (f.l == -1 || f.r == -1)
      continue;

    Vector2d flux;
    if ((flux = interpolate(f.l)).dot(f.normal) > 0)
    {
      uF(i) += flux.x();
      vF(i) += flux.y();
    }
    if ((flux = interpolate(f.r)).dot(f.normal) < 0)
    {
      uF(i) += flux.x();
      vF(i) += flux.y();
    }
  }

  const VectorXd &areas = getAreas();
  const VectorXd &nx = getNx();
  const VectorXd &ny = getNy();

  return areas.cwiseProduct(uF.cwiseProduct(nx) + vF.cwiseProduct(ny));
}

SpMat FVM::laplacian(const VectorXd &gamma)
{
  SpMat areas = convertDiagonal(getAreas());
  SpMat nx = convertDiagonal(getNx());
  SpMat ny = convertDiagonal(getNy());
}
