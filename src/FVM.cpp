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
        res(e.i) = surface.getFixed(e.i).value_or(0);
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

VectorXd FVM::calcDiv(const VectorXd &u, const VectorXd &v, ISurface &uSf, ISurface &vSf)
{
  VectorXd uF = linear(u, uSf);
  VectorXd vF = linear(v, vSf);

  VectorXd res(n);
  for (int i = 0; i < n; i++)
  {
    for (const Edge &e : grid.getAdj(i))
    {
      res(i) += Vector(uF(e.i), vF(e.i)).dot(e.n);
    }
  }
  return res;
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
      if (f.l != -1)
      {
        flux(i) = gamma(l) * f.a / surface.getDis(i, l);
      }
      if (f.r != -1)
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

        if (auto k = surface.getFixed(e.i)) {
          b(i) += *k * e.a / surface.getDis(e.i, i);
        }
      }
    }
  }

  // Debug::debugSpMat(M);

  return {M, b};
}
