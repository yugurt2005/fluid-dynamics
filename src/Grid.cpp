#include "Grid.h"

Grid::Grid(
    int n,
    int m,
    vector<Face> &faces,
    vector<double> &cellX,
    vector<double> &cellY,
    vector<std::pair<int, Vector2d>> &walls) : n(n), m(m)
{
  adjacents = new vector<int>[n];
  for (Face &face : faces)
  {
    adjacents[face.l].push_back(face.r);
    adjacents[face.r].push_back(face.l);
  }

  this->faces = faces;
  this->cellX = cellX;
  this->cellY = cellY;
  this->walls = walls;
}

Grid::~Grid()
{
  delete[] adjacents;
}

vector<int> Grid::getAdjacents(int index)
{
  return adjacents[index];
}

VectorXd Grid::getAdjDifference(int index, vector<double> &data)
{
  VectorXd answer(data.size());
  for (int i = 0; i < adjacents[index].size(); i++)
  {
    answer(i) = data[adjacents[index][i]] - data[index];
  }
  return answer;
}

VectorXd Grid::getAdjDx(int index)
{
  return getAdjDifference(index, cellX);
}

VectorXd Grid::getAdjDy(int index)
{
  return getAdjDifference(index, cellY);
}

Vector2d Grid::getCellPos(int index)
{
  return Vector2d(cellX[index], cellY[index]);
}

vector<Face> &Grid::getFaces()
{
  return faces;
}
