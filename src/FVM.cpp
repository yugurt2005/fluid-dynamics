#include "FVM.h"
#include <iostream>

FVM::FVM(int n, IGrid &grid) : n(n), grid(grid) {
  createGradientMatrix();
}

MatrixXd FVM::calcLeastSquaresGradient(int n, const VectorXd &dx, const VectorXd &dy)
{
  assert(dx.rows() == n);
  assert(dy.rows() == n);

  MatrixXd d(n, 2);
  d.col(0) = dx;
  d.col(1) = dy;

  MatrixXd w = MatrixXd::Zero(n, n);
  for (int i = 0; i < n; i++)
  {
    w(i, i) = 1 / Vector2d(dx(i), dy(i)).norm();
  }

  MatrixXd dT = d.transpose();
  MatrixXd wT = w.transpose();

  return (dT * wT * w * d).inverse() * dT * wT * w;
}

void FVM::createGradientMatrix()
{
  dxMat = SpMat(n, n);
  dyMat = SpMat(n, n);

  vector<Eigen::Triplet<double>> dxs;
  vector<Eigen::Triplet<double>> dys;

  for (int c = 0; c < n; c++)
  {
    vector<int> adjacents = getAdjacents(c);

    MatrixXd gradient = calcLeastSquaresGradient(
        adjacents.size(),
        getAdjDx(c),
        getAdjDy(c));

    for (int i = 0; i < adjacents.size(); i++)
    {
      double dx = gradient(i, 0);
      double dy = gradient(i, 1);

      dxs.emplace_back(c, adjacents[i], dx);
      dys.emplace_back(c, adjacents[i], dy);

      dxs.emplace_back(c, c, -dx);
      dys.emplace_back(c, c, -dy);
    }
  }

  dxMat.setFromTriplets(dxs.begin(), dxs.end());
  dyMat.setFromTriplets(dys.begin(), dys.end());
}
