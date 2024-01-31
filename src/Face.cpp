#include "Face.h"

Face::Face(int l, int r, Vector2d p, Vector2d q)
{
  isWall = l == -1 || r == -1;
  this->l = l;
  this->r = r;
  this->p = p;
  this->q = q;

  c = (p + q) / 2;

  Vector2d delta = p - q;

  n = Vector2d(+delta.y(), -delta.x()).normalized();

  area = delta.norm();
}