#include "FVM.h"

FVM::FVM(IGrid &grid) : grid(grid)
{
  n = grid.getN();
  z = grid.getZ();

  buildAdj();
  buildGradients();
  buildLaplacian();
}

SpMat FVM::convertDiagonal(const VectorXd &input)
{
  return input.asDiagonal().diagonal().sparseView();
}

void FVM::buildAdj()
{
  Adj = SpMat(n, z);

  const vector<Face> &faces = grid.getFaces();
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

void FVM::buildLaplacian()
{
  const vector<Face> &faces = grid.getFaces();

  laplacian = VectorXd(z);
  for (int i = 0; i < z; i++)
  {
    Face f = faces[i];

    if (f.l != -1 && f.r != -1)
    {
      laplacian(i) = std::abs(f.delta.dot(f.normal)) / f.delta.dot(f.delta);
    }
    else if (f.l != -1)
    {
      laplacian(i) = (f.center - grid.getCenter(f.l)).norm();
    }
    else if (f.r != -1)
    {
      laplacian(i) = (f.center - grid.getCenter(f.r)).norm();
    }
  }

  laplacian = laplacian.cwiseProduct(getAreas());
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

  const vector<Face> &faces = grid.getFaces();

  assert(faces.size() == z);

  VectorXd du = Gx * u;
  VectorXd dv = Gy * v;

  VectorXd uF(z);
  VectorXd vF(z);

  // Upwind Differencing
  auto interpolate = [&](int i, const Face &f)
  {
    Vector2d center = grid.getCenter(i);
    double nu = u(i) + du(i) * (f.center.x() - center.x());
    double nv = v(i) + dv(i) * (f.center.y() - center.y());
    return Vector2d(nu, nv);
  };

  for (int i = 0; i < z; i++)
  {
    Face f = faces[i];

    if (f.l == -1 || f.r == -1)
      continue;

    Vector2d flux;
    if ((flux = interpolate(f.l, f)).dot(f.normal) > 0)
    {
      uF(i) += flux.x();
      vF(i) += flux.y();
    }
    if ((flux = interpolate(f.r, f)).dot(f.normal) < 0)
    {
      uF(i) += flux.x();
      vF(i) += flux.y();
    }
  }

  const VectorXd &areas = getAreas();
  const VectorXd &nx = getNx();
  const VectorXd &ny = getNy();

  assert(uF.rows() == z);
  assert(vF.rows() == z);

  return areas.cwiseProduct(uF.cwiseProduct(nx) + vF.cwiseProduct(ny));
}

SpMat FVM::calcInterpolate(const VectorXd &flux) {
  assert(flux.size() == z);

  const vector<Face> &faces = grid.getFaces();

  SpMat dxInterpolate(z, n); // multiplies by Gx * u
  SpMat dyInterpolate(z, n); // multiplies by Gy * u
  SpMat valInterpolate(z, n); // multiplies by u

  for (int i = 0; i < z; i++) {
    Face f = faces[i];

    int node = 0;
    if (flux(i) == 0) {
      continue;
    }
    else if (flux(i) > 0) 
      node = f.l; 
    else if (flux(i) < 0) 
      node = f.r;

    if (node == -1) {
      continue;
    }

    Vector2d displacement = f.center - grid.getCenter(node);
    
    valInterpolate.insert(i, node) = 1;
    dxInterpolate.insert(i, node) = displacement.x();
    dyInterpolate.insert(i, node) = displacement.y();
  }
  SpMat xv = dxInterpolate * Gx;
  SpMat yv = dyInterpolate * Gy;
  SpMat result = valInterpolate + dxInterpolate * Gx + dyInterpolate * Gy;

  for (int i = 0; i < z; i++) {
    for (int j = 0; j < n; j++) {
      double value = result.coeff(i, j);
      std::cout << (std::abs(value) < 0.001 ? 0.0 : value) << " ";
    }
    std::cout << std::endl;
  }

  return result;
}

SpMat FVM::calcDiv(const VectorXd &fluxes)
{
  assert(fluxes.size() == z);

  return Adj * convertDiagonal(fluxes);
}

SpMat FVM::calcLaplacian(const VectorXd &gamma)
{
  assert(gamma.size() == n);

  const vector<Face> &faces = grid.getFaces();

  SpMat result(n, n);
  for (int i = 0; i < z; i++)
  {
    Face f = faces[i];

    if (f.l != -1 && f.r != -1)
    {
      double dLF = (f.center - grid.getCenter(f.l)).norm();
      double dRF = (f.center - grid.getCenter(f.r)).norm();
      double coeff = (dLF * gamma(f.l) + dRF * gamma(f.r)) / f.delta.norm() * laplacian(i);

      result.coeffRef(f.l, f.l) -= coeff;
      result.coeffRef(f.l, f.r) += coeff;

      result.coeffRef(f.r, f.r) -= coeff;
      result.coeffRef(f.r, f.l) += coeff;

      // std::cout << "Flux " << i << ": " << coeff << std::endl;
    }
    else if (f.l != -1)
    {
      result.coeffRef(f.l, f.l) -= gamma(f.l) * laplacian(i);
      // std::cout << "Flux " << i << ": " << gamma(f.l) * laplacian(i) << std::endl;
    }
    else if (f.r != -1)
    {
      result.coeffRef(f.r, f.r) -= gamma(f.r) * laplacian(i);
      // std::cout << "Flux " << i << ": " << gamma(f.r) * laplacian(i) << std::endl;
    }
    else {
      continue;
    }
  }

  for (int i = 0; i < z; i++) {
    // std::cout << i << ": " << laplacian(i) << std::endl;
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      double value = result.coeff(i, j);
      // std::cout << (std::abs(value) < 0.001 ? 0.0 : value) << " ";
    }
    // std::cout << std::endl;
  }

  return result;
}
