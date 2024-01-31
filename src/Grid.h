#ifndef GRID_H
#define GRID_H

#include <iostream>
#include <set>
#include <vector>

#include "Face.h"
#include "Edge.h"

using std::cout;
using std::endl;
using std::vector;

class Grid {
  int n, m;

  vector<Face> faces;

  vector<Edge> *adj;

  vector<Vector2d> centers;

  vector<int> boundaryLayer;

  vector<bool> isBoundaryLayer;

public:
  Grid(vector<Face> faces);

  inline int getN() const { return n; }

  inline int getM() const { return m; }

  inline const vector<Face> &getFaces() const { return faces; }

  inline const vector<Edge> &getAdj(int i) const { return adj[i]; }

  inline Vector2d getCenter(int i) const { return centers[i]; }

  Vector2d calcDis(int i, int j) const;
};

#endif // GRID_H
