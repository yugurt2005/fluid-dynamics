#include "SIMPLE.h"

SpMat calcMomentumMatrix(const State &state) {
  return {};
}

static State &relax(State &state, State &oldState, double factor) {
  state.p = state.p * factor + oldState.p * (1 - factor);
  state.k = state.k * factor + oldState.k * (1 - factor);
  state.o = state.o * factor + oldState.o * (1 - factor);
  return state;
}

State SIMPLE::step(const State &state) {
  // Calculate Velocities
  SpMat M = calcMomentumMatrix(state);

  const SpMat &Gx = getGradientMatrixX();
  const SpMat &Gy = getGradientMatrixY();

  Eigen::BiCGSTAB<SpMat> solver;
  solver.compute(M);

  VectorXd px = Gx * state.p;
  VectorXd py = Gx * state.p;

  VectorXd newU = solver.solve(px);
  VectorXd newV = solver.solve(py);

  // Define Terms
  const Eigen::Diagonal<SpMat> A = M.diagonal();
  const Eigen::Inverse<Eigen::Diagonal<SpMat>> A_I = A.inverse();

  // dim(H) = n x 1
  auto Hx = A * newU + px;
  auto Hy = A * newV + py;

  // Correct Pressure
  solver.compute(Gx * A_I * Gx + Gy * A_I * Gy);
  VectorXd p = solver.solve(Gx * A_I * Hx + Gy * A_I * Hy);

  // Correct Velocities
  VectorXd u = A_I * Hx - A_I * Gx * p;
  VectorXd v = A_I * Hy - A_I * Gy * p;

  return State{u, v, p, state.k, state.o}.relax(state, parameters.ALPHA);
}
