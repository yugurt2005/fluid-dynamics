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

  double getDis(int index) const override;

  double getDis(int index, int cell) const override;

  inline std::optional<double> getFixed(int index) const override {
    return fixed[index];
  };

  inline const vector<Face> &getFaces() const override { return faces; }
};

#endif // SURFACE_H
