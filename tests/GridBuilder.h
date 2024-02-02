#ifndef GRIDBUILDER_H
#define GRIDBUILDER_H

#include <cassert>
#include <vector>

#include "../src/Face.h"

using std::vector;

class GridBuilder {
public:
  static vector<Face> buildRectangularGrid(int n, int m, double w, double h);

  static vector<Face> buildTriangularGrid(int n, int m, double w, double h);
};

#endif // GRIDBUILDER_H
