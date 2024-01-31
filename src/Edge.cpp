#include "Edge.h"

Edge::Edge(int adj, const Face &face): adj(adj) {
  isWall = face.isWall;
  p = face.p;
  q = face.q;
  c = face.c;
  area = face.area;

  if (adj == face.r) {
    n = face.n;
  } else {
    n = -face.n;
  }
}