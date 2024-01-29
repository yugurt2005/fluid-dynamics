#ifndef CFD_FVM_H
#define CFD_FVM_H

#include <cassert>
#include <vector>

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
  int z;

  IGrid &grid;

  SpMat Adj;
  SpMat Gx;
  SpMat Gy;
  SpMat Interpolate;

  inline const VectorXd &getAreas() { return grid.getAreas(); }

  inline const VectorXd &getNx() { return grid.getNx(); }

  inline const VectorXd &getNy() { return grid.getNy(); }

  inline const vector<int> &getNeighbors(int index) { return grid.getNeighbors(index); };

public:
  FVM(IGrid &grid);

  static SpMat convertDiagonal(const VectorXd &input);

  void buildAdj();

  void buildGradients();

  void buildInterpolate();

  MatrixXd calcLeastSquaresGradient(const vector<Vector2d> &distances);

  VectorXd calcMassFlux(const VectorXd &u, const VectorXd &v);

  SpMat laplacian(const VectorXd &gamma);

  inline const SpMat &getGx() { return Gx; };

  inline const SpMat &getGy() { return Gy; };

  inline const SpMat &getAdj() { return Adj; }
};

#endif // CFD_FVM_H