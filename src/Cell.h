#ifndef CELL_H
#define CELL_H

#include <vector>

#include "Eigen/Core"

#include "Edge.h"

using Eigen::Vector2d;

struct Cell {
  int n;
  int i;

  double volume;
  Vector2d center;

  std::vector<Edge> edges;
  
  Cell(int i);

  void add();

  void init();
};

#endif // CELL_H
