#ifndef SIMPLE_H
#define SIMPLE_H

#include "Eigen/Dense"
#include "Eigen/IterativeLinearSolvers"
#include "Eigen/Sparse"

#include "State.h"
#include "Parameters.h"
#include "FVM.h"

#include "../interfaces/IHalo.h"

using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::Vector2d;
using Eigen::VectorXd;

typedef SparseMatrix<double> SpMat;

class SIMPLE
{
private:
  Parameters parameters;
  FVM &fvm;
  IHalo &halos;

  inline const SpMat &getGradientMatrixX() { return fvm.getGradientX(); }
  inline const SpMat &getGradientMatrixY() { return fvm.getGradientY(); };

  inline const vector<std::pair<int, Vector2d>> &getWalls()
  {
    return halos.getWalls();
  };

public:
  SIMPLE(Parameters parameters, FVM &fvm, IHalo &grid);

  void applyWalls(VectorXd &u, VectorXd &v);

  SpMat calculateReynoldsStressX(const State &state);
  SpMat calculateReynoldsStressY(const State &state);

  State step(const State &state);

  void loop();
};

#endif // SIMPLE_H