#include "Face.h"

Face::Face(int l, int r, Vector p, Vector q)
{
  this->l = l;
  this->r = r;

  this->p = p;
  this->q = q;

  c = (p + q) / 2;

  Vector d = p - q;
  n = Vector(d.y(), -d.x()).normalized();

  a = d.norm();
}

bool Face::isWall() const
{
  return l == -1 || r == -1;
}