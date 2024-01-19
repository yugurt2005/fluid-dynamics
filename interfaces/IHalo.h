#ifndef IHALO_H
#define IHALO_H

#include <vector>

#include "Eigen/Core"

using std::vector;

class IHalo {
public:
  virtual vector<std::pair<int, Eigen::Vector2d>> &getWalls() = 0;
};

#endif // IHALO_H
