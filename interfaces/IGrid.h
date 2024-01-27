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

  virtual int getZ() = 0;

  virtual Eigen::Vector2d getCenter(int index) = 0;

  virtual const vector<int> &getNeighbors(int index) = 0;

  virtual const Eigen::VectorXd &getAreas() = 0;

  virtual const Eigen::VectorXd &getNx() = 0;

  virtual const Eigen::VectorXd &getNy() = 0;

  virtual const vector<Face> &getFaces() = 0;
};

#endif // IGRID_H
