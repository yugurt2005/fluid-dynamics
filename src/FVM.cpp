#include "FVM.h"

FVM::FVM(const Grid &grid) : grid(grid)
{
  n = grid.getN();
  m = grid.getM();
}

VectorXd FVM::linear(const VectorXd &phi)
{
  VectorXd res(m);
  for (int i = 0; i < n; i++)
  {
    for (const Edge &e : grid.getAdj(i))
    {
      if (e.isWall)
      {
        res(e.index) = phi(i);
      }
      else
      {
        res(e.index) += phi(i) * e.d / e.dis;
      }
    }
  }
  return res;
}

std::tuple<VectorXd, VectorXd> FVM::calcDf(const VectorXd &phi)
{
  VectorXd sigma = linear(phi);

  VectorXd dx(n);
  VectorXd dy(n);

  for (int i = 0; i < n; i++)
  {
    for (const Edge &e : grid.getAdj(i))
    {
      if (!e.isWall)
      {
        dx(i) += sigma(e.index) * e.a * e.n.x();
        dy(i) += sigma(e.index) * e.a * e.n.y();
      }
      else
      {
        dx(i) = 0;
        dy(i) = 0;
      }
    }
    dx(i) /= grid.getVolume(i);
    dy(i) /= grid.getVolume(i);
  }

  return {dx, dy};
}

SpMat FVM::laplacian(const VectorXd &gamma)
{
  for (int i = 0; i < n; i++)
    assert(std::abs(gamma(i)) > 1e-9 && "The gamma coefficients must be non-zero");

  VectorXd flux(m);

  const vector<Face> &faces = grid.getFaces();
  for (int i = 0; i < m; i++)
  {
    const Face &f = faces[i];
    int l = f.l;
    int r = f.r;

    if (!f.isWall)
    {
      flux(i) = 1.0 / (f.lDel / gamma(l) + f.rDel / gamma(r)) * f.area;
    }
    else if (f.l != -1)
    {
      flux(i) = gamma(l) * f.area / f.lDel;
    }
    else if (f.r != -1)
    {
      flux(i) = gamma(r) * f.area / f.rDel;
    }
  }

  SpMat M(n, n);
  for (int i = 0; i < n; i++)
  {
    for (const Edge &e : grid.getAdj(i))
    {
      if (!e.isWall)
      {
        int a = e.to;
        M.coeffRef(i, a) += flux(e.index);
        M.coeffRef(i, i) -= flux(e.index);
      }
      else {
        M.coeffRef(i, i) -= flux(e.index);
      }
    }
  }

  // Debug::debugSpMat(M);

  return M;
}
