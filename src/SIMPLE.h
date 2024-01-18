#ifndef SIMPLE_H
#define SIMPLE_H

#include "Eigen/Dense"
#include "Eigen/IterativeLinearSolvers"
#include "Eigen/Sparse"

#include "State.h"
#include "Parameters.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SparseMatrix;
using Eigen::Triplet;

typedef SparseMatrix<double> SpMat;

class SIMPLE {
private:
  Parameters parameters;

  double getDensity();

  const SpMat &getGradientMatrixX();
  const SpMat &getGradientMatrixY();

public:
  SIMPLE(Parameters parameters);

  SpMat calculateReynoldsStress(const State &state);

  State step(const State &state);

  void loop();
};

#endif // SIMPLE_H