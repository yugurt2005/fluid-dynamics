#ifndef GRIDBUILDER_H
#define GRIDBUILDER_H

#include <cassert>

#include "../src/Grid.h"

class GridBuilder {
public:
  static Grid buildRectangularGrid(int n, int m, double w, double h);

  static Grid buildTriangularGrid(int n, int m, double w, double h);
};

#endif // GRIDBUILDER_H
