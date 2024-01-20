#ifndef IGRID_H
#define IGRID_H

#include <vector>

#include "Eigen/Core"
#include "Eigen/Dense"

#include "../src/Face.h"

using std::vector;

class IGrid
{
public:
  virtual int getN() = 0;

  virtual int getM() = 0;

  virtual vector<int> getAdjacents(int index) = 0;

  virtual Eigen::VectorXd getAdjDx(int index) = 0;

  virtual Eigen::VectorXd getAdjDy(int index) = 0;

  virtual double getCellX(int index) = 0;

  virtual double getCellY(int index) = 0;

  virtual Eigen::Vector2d getCellPos(int index) = 0;

  virtual vector<Face> &getFaces() = 0;
};

#endif // IGRID_H
