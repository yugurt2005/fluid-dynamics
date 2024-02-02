#ifndef FACE_H
#define FACE_H

#include "Eigen/Core"

typedef Eigen::Vector2d Vector;

struct Face {
  int l;
  int r;

  Vector p;
  Vector q;
  Vector c;
  Vector n;
  double a;
};

#endif // FACE_H
