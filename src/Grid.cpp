#include "Grid.h"

Grid::Grid(vector<Face> faces)
{
  this->faces = faces;

  m = faces.size();
  n = 0;
  for (Face &f : faces)
  {
    n = std::max(n, std::max(f.l, f.r) + 1);
  }

  adj = new vector<Edge>[n]();
  for (Face &f : faces)
  {
    if (f.l != -1)
    {
      adj[f.l].push_back(Edge(f.r, f));
    }
    if (f.r != -1)
    {
      adj[f.r].push_back(Edge(f.l, f));
    }
  }

  std::set<int> s;
  for (Face &f : faces)
  {
    if (f.isWall)
    {
      if (f.l != -1)
      {
        s.insert(f.l);
      }
      if (f.r != -1)
      {
        s.insert(f.r);
      }
    }
  }
  boundaryLayer = vector<int>(s.begin(), s.end());

  isBoundaryLayer = vector<bool>(n);
  for (int i : boundaryLayer)
  {
    isBoundaryLayer[i] = true;
  }

  centers = vector<Vector2d>(n);
  for (int i = 0; i < n; i++) {
    Vector2d c(0, 0);
    for (Edge &e : adj[i]) {
      c += e.c;
    }
    centers[i] = c / adj[i].size();
  }
}