#include "SIMPLE.h"

SpMat SIMPLE::calculateReynoldsStress(const State &state) {
}

State SIMPLE::step(const State &state)
{
  // Calculate Momentum Coefficients
  double alpha = parameters.ALPHA;
  double rho = parameters.RHO;
  double mu = parameters.MU;

  const auto &oldU = state.u.asDiagonal();
  const auto &oldV = state.v.asDiagonal();

  const SpMat &Gx = getGradientMatrixX();
  const SpMat &Gy = getGradientMatrixY();

  SpMat M = rho * (oldU * Gx - mu * Gx * Gx + oldV * Gy - mu * Gy * Gy);

  assert(M.rows() == M.cols() && "M is not square");

  // Compute Velocities
  Eigen::BiCGSTAB<SpMat> solver;
  solver.compute(M);

  VectorXd px = Gx * state.p;
  VectorXd py = Gx * state.p;

  VectorXd newU = solver.solve(-px);
  VectorXd newV = solver.solve(-py);

  // Define Terms
  auto A = M.diagonal();
  auto A_I = A.inverse();

  // dim(H) = [n x 1]
  auto Hx = A * newU + px;
  auto Hy = A * newV + py;

  // Correct Pressure
  solver.compute(Gx * A_I * Gx + Gy * A_I * Gy);
  VectorXd p = solver.solve(Gx * A_I * Hx + Gy * A_I * Hy);

  // Correct Velocities
  VectorXd u = A_I * Hx - A_I * Gx * p;
  VectorXd v = A_I * Hy - A_I * Gy * p;

  return State{u, v, p, state.k, state.o}.relax(state, alpha);
}
