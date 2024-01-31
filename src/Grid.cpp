#include "Grid.h"

Grid::Grid(vector<Face> faces)
{
  m = faces.size();
  n = 0;
  for (const Face &f : faces)
  {
    n = std::max(n, std::max(f.l, f.r) + 1);
  }

  std::set<int> s;
  for (const Face &f : faces)
  {
    if (f.isWall)
    {
      if (f.l != -1)
        s.insert(f.l);
      if (f.r != -1)
        s.insert(f.r);
    }
  }
  boundaryLayer = vector<int>(s.begin(), s.end());

  inBoundaryLayer = vector<bool>(n);
  for (int i : boundaryLayer)
  {
    inBoundaryLayer[i] = true;
  }

  vector<int> degree(n);

  centers = vector<Vector2d>(n);
  for (int i = 0; i < m; i++)
  {
    const Face &face = faces[i];
    int l = face.l;
    int r = face.r;

    if (l != -1)
      centers[l] += face.c, degree[l]++;
    if (r != -1)
      centers[r] += face.c, degree[r]++;
  }
  for (int i = 0; i < n; i++)
  {
    centers[i] /= degree[i];
  }

  for (int i = 0; i < m; i++)
  {
    Face &face = faces[i];
    int l = face.l;
    int r = face.r;

    if (l != -1 && r != -1)
      face.dis = (centers[l] - centers[r]).norm();
    if (l != -1)
      face.lDel = (centers[l] - face.c).norm();
    if (r != -1)
      face.rDel = (centers[r] - face.c).norm();
  }

  this->faces = faces;

  adj = new vector<Edge>[n]();
  for (int i = 0; i < m; i++)
  {
    const Face &f = faces[i];
    if (f.l != -1)
      adj[f.l].push_back(Edge(f.r, i, f));
    if (f.r != -1)
      adj[f.r].push_back(Edge(f.l, i, f));
  }

  volumes = vector<double>(n);
  for (int i = 0 ; i < n; i++)
  {
    for (const Edge &e : adj[i])
    {
      volumes[i] += e.a * e.c.dot(e.n) / 2;
    }
  }
}