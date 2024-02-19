#ifndef EDGE_H
#define EDGE_H

#include "Face.h"

struct Edge {
  int to;
  int i;

  bool side;

  Vector p;
  Vector q;
  Vector c;
  Vector n;
  double a;

  Edge();

  Edge(int to, int i, const Face &f);

  bool isWall() const;
};

#endif // EDGE_H
