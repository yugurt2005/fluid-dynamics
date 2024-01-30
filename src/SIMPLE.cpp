#include "SIMPLE.h"

#include <iostream>
#include <iomanip>
#include <cmath>

SIMPLE::SIMPLE(Parameters parameters, FVM &fvm)
    : parameters(parameters), fvm(fvm)
{
  n = fvm.getN();
  z = fvm.getZ();
}

void SIMPLE::calculateVelocity(const State &state)
{
  // TODO: volume

  double alpha = parameters.ALPHA;
  double rho = parameters.RHO;
  double mu = parameters.MU;

  const auto &oldU = state.u;
  const auto &oldV = state.v;

  const SpMat &Gx = getDxMat();
  const SpMat &Gy = getDyMat();

  VectorXd muVec = VectorXd::Constant(n, parameters.MU);

  auto fluxes = fvm.calcMassFlux(oldU, oldV);
  auto laplacian = fvm.calcLaplacian(muVec);

  assert(fluxes.rows() == z);

  auto divUU = fvm.calcDiv(fluxes);

  assert(divUU.rows() == n);
  assert(divUU.cols() == n);
  assert(laplacian.rows() == n);
  assert(laplacian.cols() == n);

  M = rho * divUU - laplacian;

  assert(M.rows() == M.cols());

  Eigen::BiCGSTAB<SpMat> solver;
  solver.compute(M);

  newU = solver.solve(-Gx * state.p);
  newV = solver.solve(-Gy * state.p);

  std::cout << "SIMPLE: Correct Velocity\n";
  VectorXd temp = Gx * state.p;
  for (int i = 0; i < n; i++) {
    std::cout << "| ";
    for (int j = 0; j < n; j++) {
      
      std::cout << std::setprecision(2) << M.coeff(i, j) << " ";
    }
    std::cout << "| * | " << newU(i) << " | = | " << std::setprecision(2) << temp(i) << " |\n";
  }
  std::cout << std::endl;
}

VectorXd SIMPLE::correctPressure(const State &state) {
  const SpMat &Gx = getDxMat();
  const SpMat &Gy = getDyMat();

  // Define Terms
  auto A = M.diagonal().asDiagonal();
  auto A_I = A.inverse();

  // dim(H) = [n x 1]
  auto Hx = A * newU - M * newU;
  auto Hy = A * newV - M * newV;

  SpMat lhs = fvm.calcLaplacian(M.diagonal().cwiseInverse());
  VectorXd temp = fvm.getAdj() * fvm.calcMassFlux(A_I * Hx, A_I * Hy);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (lhs.coeff(i, j) && std::abs(lhs.coeff(i, j)) < 1e-10) {
        lhs.coeffRef(i, j) = 0;
      }
    }
  }

  Eigen::BiCGSTAB<SpMat> solver;
  solver.compute(lhs);
  VectorXd p =
      solver.solve(temp);

  std::cout << "SIMPLE: Correct Pressure\n";
  for (int i = 0; i < n; i++) {
    std::cout << "| ";
    for (int j = 0; j < n; j++) {
      std::cout << std::setprecision(5) << lhs.coeff(i, j) << " ";
    }
    std::cout << "| * | " << std::setprecision(2) << p(i) << " | = | " << std::setprecision(2) << temp(i) << " |\n";
  }
  std::cout << std::endl;


  return p;
}

State SIMPLE::step(const State &state)
{
  // Calculate Momentum Coefficients
  double alpha = parameters.ALPHA;
  double rho = parameters.RHO;
  double mu = parameters.MU;

  calculateVelocity(state);

  // Define Terms
  const SpMat &Gx = getDxMat();
  const SpMat &Gy = getDyMat();

  auto A = M.diagonal().asDiagonal();
  auto A_I = A.inverse();

  // dim(H) = [n x 1]
  auto Hx = A * newU - M * newU;
  auto Hy = A * newV - M * newV;

  auto p = correctPressure(state);

  std::cout << "New Pressures: \n";
  for (int i = 0; i < n; i++) {
    std::cout << p(i) << " ";
  }
  std::cout << std::endl;
  std::cout << std::endl;

  // Correct Velocities
  VectorXd u = A_I * Hx - A_I * Gx * p;
  VectorXd v = A_I * Hy - A_I * Gy * p;

  std::cout << "New Velocities: \n";
  for (int i = 0; i < n; i++) {
    std::cout << u(i) << " ";
  }
  std::cout << std::endl;
  for (int i = 0; i < n; i++) {
    std::cout << v(i) << " ";
  }
  std::cout << std::endl;
  std::cout << std::endl;

  // TODO: temp, assuming 2x3
  u(0) = 1;
  u(3) = 1;

  return State{u, v, p, state.k, state.o}.relax(state, alpha);
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