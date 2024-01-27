#ifndef CFD_FVM_H
#define CFD_FVM_H

#include <cassert>
#include <vector>

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/Sparse"

#include "../interfaces/IGrid.h"

using Eigen::Vector2d;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using std::vector;

typedef SparseMatrix<double> SpMat;

class FVM
{
private:
  int n;
  int z;

  IGrid &grid;

  SpMat adj;
  SpMat Gx;
  SpMat Gy;

  inline const VectorXd &getAreas() { return grid.getAreas(); }

  inline const VectorXd &getNx() { return grid.getNx(); }

  inline const VectorXd &getNy() { return grid.getNy(); }

  inline const vector<int> &getNeighbors(int index);

public:
  FVM(IGrid &grid);

  static SpMat convertDiagonal(const VectorXd &input);

  void buildAdj();

  void buildGradients();

  MatrixXd calcLeastSquaresGradient(const vector<Vector2d> &distances);

  VectorXd calcMassFlux(const VectorXd &u, const VectorXd &v);

  SpMat calcInterpolate();

  SpMat div(const VectorXd &flux);

  SpMat laplacian();

  inline const SpMat &getGx() { return Gx; };

  inline const SpMat &getGy() { return Gy; };

  inline const SpMat &getAdj() { return adj; }
};

#endif // CFD_FVM_H