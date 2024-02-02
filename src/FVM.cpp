#include "FVM.h"

FVM::FVM(IGrid &grid) : grid(grid)
{
  n = grid.getN();
  m = grid.getM();
}

VectorXd FVM::linear(const VectorXd &phi, ISurface &surface)
{
  VectorXd res(m);
  for (int i = 0; i < n; i++)
  {
    for (const Edge &e : grid.getAdj(i))
    {
      if (e.isWall())
      {
        res(e.i) = surface.getFixed(e.i).value_or(phi(i));
      }
      else
      {
        res(e.i) += phi(i) * surface.getDis(e.i, i) / surface.getDis(e.i);
      }
    }
  }
  return res;
}

std::tuple<VectorXd, VectorXd> FVM::calcDf(const VectorXd &phi, ISurface &surface)
{
  VectorXd sigma = linear(phi, surface);

  VectorXd dx(n);
  VectorXd dy(n);

  for (int i = 0; i < n; i++)
  {
    for (const Edge &e : grid.getAdj(i))
    {
      dx(i) += sigma(e.i) * e.a * e.n.x();
      dy(i) += sigma(e.i) * e.a * e.n.y();
    }
    dx(i) /= grid.getVolume(i);
    dy(i) /= grid.getVolume(i);
  }

  return {dx, dy};
}

std::tuple<SpMat, VectorXd> FVM::laplacian(const VectorXd &gamma, ISurface &surface)
{
  for (int i = 0; i < n; i++)
    assert(std::abs(gamma(i)) > 1e-9 && "The gamma coefficients must be non-zero");

  SpMat M(n, n);
  VectorXd b(n);

  VectorXd flux(m);

  const std::vector<Face> &faces = surface.getFaces();
  for (int i = 0; i < m; i++)
  {
    const Face &f = faces[i];
    int l = f.l;
    int r = f.r;

    if (!f.isWall())
    {
      double ld = surface.getDis(i, l);
      double rd = surface.getDis(i, r);
      flux(i) = 1.0 / (ld / gamma(l) + rd / gamma(r)) * f.a;
    }
    else
    {
      if (!surface.getFixed(i).has_value())
      {
        flux(i) = 0;
      }
      else if (f.l != -1)
      {
        flux(i) = gamma(l) * f.a / surface.getDis(i, l);
      }
      else if (f.r != -1)
      {
        flux(i) = gamma(r) * f.a / surface.getDis(i, r);
      }
    }
  }

  for (int i = 0; i < n; i++)
  {
    for (const Edge &e : grid.getAdj(i))
    {
      int a = e.to;
      if (!e.isWall())
      {
        M.coeffRef(i, a) += flux(e.i);
        M.coeffRef(i, i) -= flux(e.i);
      }
      else
      {
        M.coeffRef(i, i) -= flux(e.i);

        double k = surface.getFixed(e.i).value_or(0);
        if (k) {
          b(i) -= k * e.a / surface.getDis(e.i, i);
        }
      }
    }
  }

  // Debug::debugSpMat(M);

  return {M, b};
}
