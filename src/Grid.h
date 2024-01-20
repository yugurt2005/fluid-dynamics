#ifndef GRID_H
#define GRID_H

#include "../interfaces/IGrid.h"
#include "../interfaces/IHalo.h"

using Eigen::Vector2d;
using Eigen::VectorXd;

class Grid : public IGrid, public IHalo
{
private:
  int n;
  int m;

  vector<int> *adjacents;

  vector<Face> faces;

  vector<double> cellX;
  vector<double> cellY;

  vector<std::pair<int, Vector2d>> walls;

  VectorXd getAdjDifference(int index, vector<double> &data);

public:
  Grid(
      int n,
      int m,
      vector<Face> &faces,
      vector<double> &cellX,
      vector<double> &cellY,
      vector<std::pair<int, Vector2d>> &walls);

  ~Grid();

  inline int getN() override { return n; }

  inline int getM() override { return m; }

  vector<int> getAdjacents(int index) override;

  VectorXd getAdjDx(int index) override;

  VectorXd getAdjDy(int index) override;

  inline double getCellX(int index) override { return cellX[index]; }

  inline double getCellY(int index) override { return cellY[index]; }

  Vector2d getCellPos(int index) override;

  vector<Face> &getFaces() override;

  vector<std::pair<int, Eigen::Vector2d>> &getWalls() override;
};

#endif // GRID_H