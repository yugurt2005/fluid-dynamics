#include "FVM.h"

FVM::FVM(IGrid &grid) : grid(grid)
{
  n = grid.getN();
  m = grid.getM();
}

void FVM::init()
{
  createGradientMatrix();
}

VectorXd FVM::calcLeastSquaresGradient(int n, const VectorXd &dx, const VectorXd &dy)
{
  MatrixXd d(n, 2);
  d.col(0) = dx;
  d.col(1) = dy;

  MatrixXd w =
      (dx.array().square() * dy.array().square()).rsqrt().matrix().asDiagonal();

  MatrixXd dT = d.transpose();
  MatrixXd wT = w.transpose();

  return (dT * wT * w * d).inverse() * dT * wT * w;
}

void FVM::createGradientMatrix()
{
  gradientX = SpMat(n, n);
  gradientY = SpMat(n, n);

  vector<Eigen::Triplet<double>> xTrips;
  vector<Eigen::Triplet<double>> yTrips;

  for (int c = 0; c < n; c++)
  {
    vector<int> adjacents = getAdjacents(c);

    VectorXd gradient = calcLeastSquaresGradient(
        adjacents.size(),
        getAdjDx(c),
        getAdjDy(c));

    for (int i = 0; i < adjacents.size(); i++)
    {
      double dx = gradient(i);
      double dy = gradient(i);

      xTrips.emplace_back(c, adjacents[i], dx);
      yTrips.emplace_back(c, adjacents[i], dy);

      xTrips.emplace_back(c, c, -dx);
      yTrips.emplace_back(c, c, -dy);
    }
  }

  gradientX.setFromTriplets(xTrips.begin(), xTrips.end());
  gradientY.setFromTriplets(yTrips.begin(), yTrips.end());
}

double FVM::interpolateAt(int index, Face &face, double phi, double dx, double dy)
{
  return phi + (face.center - grid.getCellPos(index)).dot(Vector2d(dx, dy));
}

VectorXd FVM::interpolate(VectorXd &phi, VectorXd &flux)
{
  VectorXd dx = getGradientX() * phi;
  VectorXd dy = getGradientY() * phi;

  VectorXd answer = VectorXd::Zero(m);

  vector<Face> &faces = getFaces();
  for (int i = 0; i < m; i++)
  {
    Face &face = faces[i];

    int l = face.l;
    int r = face.r;

    if (flux(i) > 0)
      answer(i) = interpolateAt(l, face, phi(l), dx(l), dy(l));
    if (flux(i) < 0)
      answer(i) = interpolateAt(r, face, phi(r), dx(r), dy(r));
  }

  return answer;
}

VectorXd FVM::calcMassFlux(const VectorXd &us, const VectorXd &vs)
{
  VectorXd phi = VectorXd::Zero(m);

  VectorXd dxU = getGradientX() * us;
  VectorXd dyU = getGradientY() * us;
  VectorXd dxV = getGradientX() * vs;
  VectorXd dyV = getGradientY() * vs;

  auto calcMassFluxAt = [&](int index, Face &face, double u, double v) {
    double nu = interpolateAt(index, face, u, dxU(index), dyU(index));
    double nv = interpolateAt(index, face, v, dxV(index), dyV(index));

    return Vector2d(nu, nv).dot(face.normal) * face.length;
  };

  vector<Face> &faces = getFaces();
  for (int i = 0; i < m; i++)
  {
    Face &face = faces[i];

    int l = face.l;
    int r = face.r;

    if (l != -1 && Vector2d(us(l), vs(l)).dot(face.normal) > 0)
    {
      phi(i) += calcMassFluxAt(l, face, us(l), vs(l));
    }
    if (r != -1 && Vector2d(us(r), vs(r)).dot(face.normal) < 0)
    {
      phi(i) += calcMassFluxAt(r, face, us(r), vs(r));
    }
  }

  return phi;
}