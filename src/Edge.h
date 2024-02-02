#ifndef EDGE_H
#define EDGE_H

#include "Face.h"

struct Edge {
  int to;
  int i;

  Vector p;
  Vector q;
  Vector c;
  Vector n;
  double a;

  Edge();

  Edge(int to, int i, const Face &f);
};

#endif // EDGE_H
