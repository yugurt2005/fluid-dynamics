#include "Grid.h"

Grid::Grid(vector<Vector2d> centers, vector<Face> faces)
{
  this->centers = centers;

  n = centers.size();

  neighbors = new vector<int>[n]();
  for (Face face : faces)
  {
    int u = face.l;
    int v = face.r;

    assert(u < n);
    assert(v < n);

    if (u > -1 && v > -1)
    {
      neighbors[u].push_back(v);
      neighbors[v].push_back(u);
    }
  }

  z = faces.size();

  areas = VectorXd(z);
  nx = VectorXd(z);
  ny = VectorXd(z);

  for (int i = 0; i < z; i++)
  {
    Face face = faces[i];

    areas(i) = face.area;
    nx(i) = face.normal.x();
    ny(i) = face.normal.y();
  }

  for (int i = 0; i < z; i++) {
    Face &f = faces[i];
    if (f.l == -1 || f.r == -1)
      continue;
    f.updateDelta(centers[f.r] - centers[f.l]);
  }

  this->faces = faces;
}

Grid::~Grid()
{
  delete[] neighbors;
}