#ifndef SIMPLE_H
#define SIMPLE_H

#include <cassert>

#include "Eigen/Dense"
#include "Eigen/IterativeLinearSolvers"
#include "Eigen/Sparse"

#include "State.h"
#include "Parameters.h"
#include "FVM.h"

#include "../interfaces/IHalo.h"

using Eigen::MatrixXd;
using Eigen::Vector2d;
using Eigen::VectorXd;
using Eigen::SparseMatrix;

typedef SparseMatrix<double> SpMat;

class SIMPLE
{
private:
  Parameters parameters;

  int n;

  FVM &fvm;
  IHalo &grid;

  VectorXd newU;
  VectorXd newV;
  SpMat M;

  const SpMat &getDxMat() { return fvm.getDxMat(); }

  const SpMat &getDyMat() { return fvm.getDyMat(); }

  inline const vector<std::pair<int, Vector2d>> &getWalls() {
    return grid.getWalls();
  };

public:
  SIMPLE(Parameters parameters, FVM &fvm, IHalo &grid);

  inline VectorXd &getNewU() { return newU; }

  inline VectorXd &getNewV() { return newV; }

  void applyWalls(VectorXd &u, VectorXd &v);

  SpMat calculateReynoldsStressX(const State &state);

  SpMat calculateReynoldsStressY(const State &state);

  void calculateVelocity(const State &state);

  VectorXd correctPressure(const State &state);

  State step(const State &state);

  State loop(const State &state, int iterations);
};

#endif // SIMPLE_H