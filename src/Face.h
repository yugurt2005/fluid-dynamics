#ifndef FACE_H
#define FACE_H

#include <Eigen/Core>

using Eigen::Vector2d;

/*
  p
l | r
  q
*/
struct Face {
  bool isWall;
  int l;
  int r;
  Vector2d p;
  Vector2d q;
  Vector2d c;
  Vector2d n;
  double area;

  double dis;
  double lDel;
  double rDel;

  Face(int l, int r, Vector2d p, Vector2d q);
};

#endif // FACE_H
