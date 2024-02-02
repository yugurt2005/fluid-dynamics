#ifndef FVM_H
#define FVM_H

#include <cassert>
#include <tuple>

#include "Eigen/Dense"
#include "Eigen/Sparse"

#include "../interfaces/IGrid.h"
#include "../interfaces/ISurface.h"

using Eigen::VectorXd;

typedef Eigen::SparseMatrix<double> SpMat;

class FVM {
  int n;
  int m;
  IGrid &grid;

public:
  FVM(IGrid &grid);

  VectorXd linear(const VectorXd &phi, ISurface &surface);

  std::tuple<VectorXd, VectorXd> calcDf(const VectorXd &phi, ISurface &surface);

  std::tuple<SpMat, VectorXd> laplacian(const VectorXd &gamma, ISurface &surface);
};

#endif // FVM_H
