#ifndef FACE_H
#define FACE_H

#include "Eigen/Core"

using Eigen::Vector2d;

struct Face
{
  Vector2d a;
  Vector2d b;
  Vector2d center;
  Vector2d normal;
  double length;
};

#endif // FACE_H
