#ifndef IGRID_H
#define IGRID_H

#include <vector>

#include "Eigen/Core"

#include "../src/Edge.h"

typedef Eigen::Vector2d Vector;

class IGrid {
public:
  virtual int getN() = 0;

  virtual int getM() = 0;

  virtual double getVolume(int index) = 0;

  virtual Vector getCenter(int index) = 0;

  virtual const std::vector<Edge> &getAdj(int index) = 0;
};

#endif // IGRID_H
