#include "Face.h"

Face::Face(int l, int r, Vector p, Vector q)
{
  this->l = l;
  this->r = r;

  this->p = p;
  this->q = q;

  c = (p + q) / 2;

  Vector d = q - p;
  n = {d.y(), -d.x()};

  a = d.norm();
}

bool Face::isWall() const
{
  return l == -1 || r == -1;
}