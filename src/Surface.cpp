#include "Surface.h"

Surface::Surface(vector<std::pair<int, double>> fixed) {
  this->fixed.resize(n, std::nullopt);
  for (auto &f : fixed) {
    this->fixed[f.first] = f.second;
  }
}

void Surface::init(vector<Face> faces, IGrid &grid) {
  n = grid.getN();
  m = grid.getN();

  dis = vector<double>(m);
  disL = vector<double>(m);
  disR = vector<double>(m);

  for (int i = 0; i < m; i++) {
    const Face &f = faces[i];
    if (f.l != -1 && f.r != -1) {
      dis[i] = (grid.getCenter(f.l) - grid.getCenter(f.r)).norm();
    }
    if (f.l != -1) {
      disL[i] = (grid.getCenter(f.l) - f.c).norm();
    }
    if (f.r != -1) {
      disR[i] = (grid.getCenter(f.r) - f.c).norm();
    }
  }

  Surface::faces = faces;
}

double Surface::getDis(int index) { return dis[index]; }

double Surface::getDis(int index, int cell) {
  if (cell == faces[index].l) {
    return disL[index];
  }
  if (cell == faces[index].r) {
    return disR[index];
  }

  assert(false && "Invalid Indices!");
}
