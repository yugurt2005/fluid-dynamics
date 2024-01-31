#include "Edge.h"

Edge::Edge(int to, int index, const Face &face)
{
  this->to = to;
  this->index = index;

  isWall = face.isWall;

  p = face.p;
  q = face.q;
  c = face.c;
  a = face.area;

  if (to == face.r)
  {
    n = face.n;
    d = face.lDel;
  }
  else
  {
    n = -face.n;
    d = face.rDel;
  }
  
  dis = face.dis;
}