#include "SIMPLE.h"

SIMPLE::SIMPLE(Parameters parameters, FVM &fvm, IHalo &grid)
    : parameters(parameters), fvm(fvm), grid(grid) {}

void SIMPLE::applyWalls(VectorXd &u, VectorXd &v)
{
  const auto &walls = getWalls();

  for (const auto &wall : walls)
  {
    u(wall.first) = wall.second.x();
    v(wall.first) = wall.second.y();
  }
}

State SIMPLE::step(const State &state)
{
  // // Calculate Momentum Coefficients
  // double alpha = parameters.ALPHA;
  // double rho = parameters.RHO;
  // double mu = parameters.MU;

  // const auto &oldU = state.u.asDiagonal();
  // const auto &oldV = state.v.asDiagonal();

  // const SpMat &Gx = getDxMat();
  // const SpMat &Gy = getDyMat();

  // SpMat M = rho * (oldU * Gx - mu * Gx * Gx + oldV * Gy - mu * Gy * Gy);

  // // Compute Velocities
  // Eigen::BiCGSTAB<SpMat> solver;
  // solver.compute(M);

  // VectorXd newU = solver.solve(Gx * state.p);
  // VectorXd newV = solver.solve(Gy * state.p);
  // applyWalls(newU, newV);

  // // Define Terms
  // auto A = M.diagonal();
  // auto A_I = A.inverse();

  // // dim(H) = [n x 1]
  // auto Hx = A * newU - M * newU;
  // auto Hy = A * newV - M * newV;

  // // Correct Pressure
  // solver.compute(Gx * A_I * Gx + Gy * A_I * Gy);
  // VectorXd p =
  //     solver.solve(Gx * A_I * Hx + Gy * A_I * Hy);

  // // Correct Velocities
  // VectorXd u = A_I * Hx - A_I * Gx * p;
  // VectorXd v = A_I * Hy - A_I * Gy * p;

  // return State{u, v, p, state.k, state.o}.relax(state, alpha);

  return state;
}

State SIMPLE::loop(const State &state, int iterations)
{
  State result = state;
  while (iterations--)
  {
    result = step(result);
  }
  return result;
}