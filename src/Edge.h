#ifndef EDGE_H
#define EDGE_H

#include <Eigen/Core>

#include "Face.h"

using Eigen::Vector2d;

struct Edge {
  int to;
  int index;

  bool isWall;

  Vector2d p;
  Vector2d q;
  Vector2d c;
  Vector2d n;
  double a;

  double d;
  double dis;

  Edge(int to, int index, const Face &face);
};

#endif // EDGE_H
