#ifndef GRID_H
#define GRID_H

#include <vector>

#include "../interfaces/IGrid.h"

#include "Edge.h"
#include "Face.h"

using std::vector;

class Grid : public IGrid {
  int n, m;

  double *volume;
  Vector *center;

  vector<Edge> *edges;

public:
  Grid(const vector<Face> &faces);

  ~Grid();

  static double calcVolume(const vector<Edge> &faces);

  static Vector calcCenter(const vector<Edge> &faces, double volume);

  inline int getN() override { return n; }

  inline int getM() override { return m; }

  inline double getVolume(int index) override { return volume[index]; }

  inline Vector getCenter(int index) override { return center[index]; }

  inline const vector<Edge> &getAdj(int index) override { return edges[index]; }
};

#endif // GRID_H
