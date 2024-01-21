#ifndef GRID_H
#define GRID_H

#include "Face.h"

#include "../interfaces/IGrid.h"
#include "../interfaces/IHalo.h"

using Eigen::Vector2d;
using Eigen::VectorXd;

class Grid : public IGrid, public IHalo
{
private:
  int n;

  vector<int> *adjacents;

  vector<Face> faces;

  vector<double> cellX;
  vector<double> cellY;

  vector<std::pair<int, Vector2d>> walls;

  VectorXd getAdjDifference(int index, vector<double> &data);

public:
  Grid(
      int n,
      vector<std::pair<int, int>> &connections,
      vector<double> &cellX,
      vector<double> &cellY,
      vector<std::pair<int, Vector2d>> &walls);

  ~Grid();

  vector<int> getAdjacents(int index) override;

  VectorXd getAdjDx(int index) override;

  VectorXd getAdjDy(int index) override;

  vector<std::pair<int, Vector2d>> &getWalls() override;
};

#endif // GRID_H