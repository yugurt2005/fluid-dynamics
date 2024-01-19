#include "Grid.h"

Grid::Grid(
    int n,
    vector<std::pair<int, int>> connections,
    vector<double> &cellX,
    vector<double> &cellY,
    vector<std::pair<int, Vector2d>> &walls)
{
  adjacents = new vector<int>[n];
  for (auto [u, v] : connections)
  {
    adjacents[u].push_back(v);
    adjacents[v].push_back(u);
  }

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

VectorXd Grid::getDifference(int index, vector<double> &data)
{
  VectorXd answer(data.size());
  for (int i = 0; i < data.size(); i++)
  {
    answer(i) = data[i] - data[index];
  }
  return answer;
}

VectorXd Grid::getAdjDx(int index)
{
  return getDifference(index, cellX);
}

VectorXd Grid::getAdjDy(int index)
{
  return getDifference(index, cellY);
}
