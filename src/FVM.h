#ifndef FVM_H
#define FVM_H

#include <cassert>
#include <tuple>

#include "Eigen/Dense"
#include "Eigen/Sparse"

#include "Debug.h"
#include "../interfaces/IGrid.h"
#include "../interfaces/ISurface.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;

using std::vector;

typedef Eigen::SparseMatrix<double> SpMat;

class FVM {
  int n;
  int m;
  IGrid &grid;

public:
  FVM(IGrid &grid);

  VectorXd linear(const VectorXd &phi, ISurface &surface);

  std::tuple<SpMat, VectorXd, SpMat, VectorXd> buildGradients(ISurface &surface);

  std::tuple<VectorXd, VectorXd> calcDf(const VectorXd &phi, ISurface &surface);

  VectorXd calcDiv(const VectorXd &u, const VectorXd &v, ISurface &uSf, ISurface &vSf);

  std::tuple<SpMat, VectorXd> laplacian(const VectorXd &gamma, ISurface &surface);

  MatrixXd calcLeastSquaresGradient(const vector<Vector> &distances);
};

#endif // FVM_H
