#ifndef GRID_H
#define GRID_H

#include "Face.h"

#include "../interfaces/IGrid.h"

using Eigen::Vector2d;
using Eigen::VectorXd;

class Grid : public IGrid
{
private:
  int n;
  int z;

  vector<int> *neighbors;

  vector<Vector2d> centers;

  vector<Face> faces;
  VectorXd areas;
  VectorXd nx;
  VectorXd ny;

public:
  Grid(vector<Vector2d> centers, vector<Face> faces);
      
  ~Grid();

  inline int getN() override { return n; };

  inline int getZ() override { return z; };

  inline Vector2d getCenter(int index) override { return centers[index]; };

  inline const vector<int> &getNeighbors(int index) override { return neighbors[index]; };

  inline const VectorXd &getAreas() override { return areas; };

  inline const VectorXd &getNx() override { return nx; };

  inline const VectorXd &getNy() override { return ny; };

  inline const vector<Face> &getFaces() override { return faces; };
};

#endif // GRID_H