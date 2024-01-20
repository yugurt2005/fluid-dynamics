#include "SIMPLE.h"

SIMPLE::SIMPLE(Parameters parameters, FVM &fvm, IHalo &grid)
    : parameters(parameters), fvm(fvm), halos(grid) {}

void SIMPLE::applyWalls(VectorXd &u, VectorXd &v)
{
  const auto &walls = getWalls();

  for (const auto &wall : walls)
  {
    u(wall.first) = wall.second.x();
    v(wall.first) = wall.second.y();
  }
}

SpMat SIMPLE::calculateReynoldsStressX(const State &state) {}

SpMat SIMPLE::calculateReynoldsStressY(const State &state) {}

State SIMPLE::step(const State &state)
{
  
}
