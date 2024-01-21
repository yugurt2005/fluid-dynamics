#ifndef CFD_FVM_H
#define CFD_FVM_H

#include <vector>
#include <cassert>

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/Sparse"

#include "../interfaces/IGrid.h"

using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::Vector2d;
using Eigen::VectorXd;
using std::vector;

typedef SparseMatrix<double> SpMat;

class FVM
{
private:
  int n;
  IGrid &grid;
  SparseMatrix<double> dxMat;
  SparseMatrix<double> dyMat;

  inline vector<int> getAdjacents(int index) { return grid.getAdjacents(index); }

  inline VectorXd getAdjDx(int index) { return grid.getAdjDx(index); }

  inline VectorXd getAdjDy(int index) { return grid.getAdjDy(index); };

public:
  FVM(int n, IGrid &grid);

  static MatrixXd calcLeastSquaresGradient(int n, const VectorXd &dx, const VectorXd &dy);

  void createGradientMatrix();

  inline const SpMat &getDxMat() { return dxMat; };

  inline const SpMat &getDyMat() { return dyMat; };
};

#endif // CFD_FVM_H