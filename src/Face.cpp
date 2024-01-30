#include "Face.h"

/*
    b
| L | R |
    a
*/
Face::Face(int l, int r, Vector2d a, Vector2d b)
{
  this->l = l;
  this->r = r;
  this->a = a;
  this->b = b;

  Vector2d vector = (b - a);

  this->center = (a + b) / 2;
  this->normal = Vector2d(+vector.y(), -vector.x()).normalized();
  this->area = vector.norm();
}

void Face::updateDelta(Vector2d d) {
  delta = d;
}