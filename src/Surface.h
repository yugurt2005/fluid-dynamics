#ifndef SURFACE_H
#define SURFACE_H

#include "../interfaces/ISurface.h"

#include "Grid.h"

class Surface : public ISurface {
  static int n, m;

  static vector<Face> faces;

  static vector<double> dis;
  static vector<double> disL;
  static vector<double> disR;

  static bool initialized;

  vector<std::optional<double>> fixed;

public:
  Surface(vector<std::pair<int, double>> fixed);

  static void init(vector<Face> faces, IGrid &grid);

  double getDis(int index) override;

  double getDis(int index, int cell) override;

  inline std::optional<double> getFixed(int index) override {
    return fixed[index];
  };

  inline const vector<Face> &getFaces() override { return faces; }
};

#endif // SURFACE_H
