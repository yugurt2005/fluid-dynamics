#ifndef IGRID_H
#define IGRID_H

#include <vector>

#include "Eigen/Core"
#include "Eigen/Dense"

using std::vector;

class IGrid
{
public:
  virtual vector<int> getAdjacents(int index) = 0;

  virtual Eigen::VectorXd getAdjDx(int index) = 0;

  virtual Eigen::VectorXd getAdjDy(int index) = 0;
};

#endif // IGRID_H
