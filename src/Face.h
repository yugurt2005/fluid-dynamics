#ifndef FACE_H
#define FACE_H

#include "Eigen/Core"

using Eigen::Vector2d;

struct Face
{
  int l;
  int r;
  Vector2d a;
  Vector2d b;
  Vector2d center;
  Vector2d normal;
  Vector2d delta;
  double area;

  Face(int l, int r, Vector2d a, Vector2d b);

  void updateDelta(Vector2d d);
};

#endif // FACE_H
