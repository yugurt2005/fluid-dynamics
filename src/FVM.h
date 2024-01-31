#ifndef FVM_H
#define FVM_H

#include <cassert>
#include <tuple>

#include "Eigen/Dense"
#include "Eigen/Sparse"

#include "Grid.h"
#include "Debug.h"

using Eigen::VectorXd;

typedef Eigen::SparseMatrix<double> SpMat;

class FVM {
  int n;
  int m;
  const Grid &grid;

public:
  FVM(const Grid &grid);

  VectorXd linear(const VectorXd &phi);

  std::tuple<VectorXd, VectorXd> calcDf(const VectorXd &phi);

  SpMat laplacian(const VectorXd &gamma);
};

#endif // FVM_H
