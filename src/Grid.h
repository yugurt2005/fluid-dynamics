#ifndef GRID_H
#define GRID_H

#include <iostream>
#include <set>
#include <vector>

#include "Face.h"
#include "Edge.h"

using std::vector;

class Grid {
  int n, m;

  vector<Face> faces;

  vector<Edge> *adj;

  vector<Vector2d> centers;

  vector<int> boundaryLayer;

  vector<bool> inBoundaryLayer;

  vector<double> volumes;

public:
  Grid(vector<Face> faces);

  inline int getN() const { return n; }

  inline int getM() const { return m; }

  inline const vector<Face> &getFaces() const { return faces; }

  inline const vector<Edge> &getAdj(int i) const { return adj[i]; }

  inline Vector2d getCenter(int i) const { return centers[i]; }

  inline bool isInBoundaryLayer(int i) const { return inBoundaryLayer[i]; }

  inline bool isWall(int i) const { return faces[i].isWall; }

  inline double getVolume(int i) const { return volumes[i]; };
};

#endif // GRID_H
