#ifndef SIMPLE_H
#define SIMPLE_H

#include "Eigen/Dense"
#include "Eigen/IterativeLinearSolvers"
#include "Eigen/Sparse"
#include <tuple>

#include "State.h"
#include "Parameters.h"

using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::VectorXd;

typedef SparseMatrix<double> SpMat;

class SIMPLE {
private:
  Parameters parameters;

  double getDensity();

  const SpMat &getGradientMatrixX();
  const SpMat &getGradientMatrixY();

  SpMat calcMomentumMatrix(const State &state);

public:
  SIMPLE(Parameters parameters);

  State step(const State &state);

  void loop();
};

#endif // SIMPLE_H