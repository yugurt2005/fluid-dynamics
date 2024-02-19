#include "FVM.h"

FVM::FVM(IGrid &grid) : grid(grid)
{
  n = grid.getN();
  m = grid.getM();
}

std::tuple<SpMat, VectorXd, SpMat, VectorXd> FVM::buildGradients(ISurface &surface)
{
  SpMat Gx = SpMat(n, n);
  SpMat Gy = SpMat(n, n);

  VectorXd bx = VectorXd::Zero(n);
  VectorXd by = VectorXd::Zero(n);

  for (int c = 0; c < n; c++)
  {
    vector<Vector> distances;
    for (Edge e : grid.getAdj(c))
    {
      if (!e.isWall())
      {
        distances.push_back(grid.getCenter(e.to) - grid.getCenter(c));
      }
      else
      {
        distances.push_back(e.c - grid.getCenter(c));
      }
    }

    for (int i = 0; i < distances.size(); i++)
    {
      MatrixXd gradient = calcLeastSquaresGradient(distances);

      Edge adj = grid.getAdj(c)[i];
      if (adj.isWall())
      {
        bx(c) += gradient(0, i) * surface.getFixed(adj.i).value_or(0);
        by(c) += gradient(1, i) * surface.getFixed(adj.i).value_or(0);

        Gx.coeffRef(c, c) -= gradient(0, i);
        Gy.coeffRef(c, c) -= gradient(1, i);
      }
      else
      {
        Gx.coeffRef(c, adj.to) += gradient(0, i);
        Gy.coeffRef(c, adj.to) += gradient(1, i);

        Gx.coeffRef(c, c) -= gradient(0, i);
        Gy.coeffRef(c, c) -= gradient(1, i);
      }
    }
  }

  return {Gx, bx, Gy, by};
}

MatrixXd FVM::calcLeastSquaresGradient(const vector<Vector> &distances)
{
  int n = distances.size();

  MatrixXd d(n, 2);
  for (int i = 0; i < n; i++)
  {
    d(i, 0) = distances[i].x();
    d(i, 1) = distances[i].y();
  }

  MatrixXd w = MatrixXd::Zero(n, n);
  for (int i = 0; i < n; i++)
  {
    w(i, i) = 1.0 / distances[i].norm();
  }

  MatrixXd dT = d.transpose();
  MatrixXd wT = w.transpose();

  return (dT * wT * w * d).inverse() * dT * wT * w;
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

        if (auto k = surface.getFixed(e.i))
        {
          b(i) += *k * e.a / surface.getDis(e.i, i);
        }
      }
    }
  }

  // Debug::debugSpMat(M);

  return {M, b};
}

std::tuple<SpMat, VectorXd> FVM::calcConvection(const VectorXd &flux, ISurface &surface)
{
  SpMat M(n, n);
  VectorXd b(n);
  for (int i = 0; i < n; i++)
  {
    for (const Edge &e : grid.getAdj(i))
    {
      if (!e.isWall())
      {
        M.coeffRef(i, e.to) += flux(e.i);
      }
      else
      {
        b(i) += flux(e.i);
      }
    }
  }

  return {M, b};
}

VectorXd FVM::calcMassFlux(const VectorXd &u, const VectorXd &v, ISurface &uSf, ISurface &vSf)
{
  assert(u.size() == n);
  assert(v.size() == n);

  VectorXd flux(m);

  auto [Gux, bux, Guy, buy] = buildGradients(uSf);
  auto [Gvx, bvx, Gvy, bvy] = buildGradients(vSf);

  VectorXd dux = Gux * u + bux;
  VectorXd duy = Guy * u + buy;
  VectorXd dvx = Gvx * v + bvx;
  VectorXd dvy = Gvy * v + bvy;

  // Debug::debug2d(dux, 2, 2, "dux");
  // Debug::debug2d(duy, 2, 2, "duy");
  // Debug::debug2d(dvx, 2, 2, "dvx");
  // Debug::debug2d(dvy, 2, 2, "dvy");

  // Upwind Differencing
  auto interpolate = [&](int i, Vector c)
  {
    Vector center = grid.getCenter(i);
    double nu = u(i) + (c - center).dot(Vector(dux(i), duy(i)));
    double nv = v(i) + (c - center).dot(Vector(dvx(i), dvy(i)));
    return Vector(nu, nv);
  };

  for (int i = 0; i < n; i++)
  {
    for (const Edge &e : grid.getAdj(i))
    {
      if (e.isWall())
      {
        double u = uSf.getFixed(e.i).value_or(0);
        double v = vSf.getFixed(e.i).value_or(0);
        Vector n = e.side ? e.n : -e.n;
        flux(e.i) = Vector(u, v).dot(n);
      }
      else
      {
        Vector velocity = interpolate(i, e.c);

        double f;
        if ((f = velocity.dot(e.n)) > 0)
        {
          if (e.side)
          {
            flux(e.i) += f * e.a;
          }
          else
          {
            flux(e.i) -= f * e.a;
          }
        }
      }
    }
  }

  return flux;
}