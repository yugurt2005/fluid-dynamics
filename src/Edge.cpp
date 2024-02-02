#include "Edge.h"

Edge::Edge(int to, int i, const Face &f)
{
  this->to = to;
  this->i = i;

  p = f.p;
  q = f.q;
  c = f.c;
  a = f.a;

  if (f.r == to)
    n = f.n;
  if (f.l == to)
    n = f.n * -1;
}

bool Edge::isWall() const
{
  return to == -1;
}