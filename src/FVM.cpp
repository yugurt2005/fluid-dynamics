#include "FVM.h"

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

  vector<Eigen::Triplet<double>> dxs;
  vector<Eigen::Triplet<double>> dys;

  for (int c = 0; c < n; c++)
  {
    vector<int> adjacents = getAdjacents(c);

    VectorXd gradient = calcLeastSquaresGradient(
      adjacents.size(),
      getAdjDx(c),
      getAdjDy(c)
    );

    for (int i = 0; i < adjacents.size(); i++)
    {
      double dx = gradient(i);
      double dy = gradient(i);

      dxs.emplace_back(c, adjacents[i], dx);
      dys.emplace_back(c, adjacents[i], dy);

      dxs.emplace_back(c, c, -dx);
      dys.emplace_back(c, c, -dy);
    }
  }

  gradientX.setFromTriplets(dxs.begin(), dxs.end());
  gradientY.setFromTriplets(dys.begin(), dys.end());
}
