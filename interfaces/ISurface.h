#ifndef ISURFACE_H
#define ISURFACE_H

#include <optional>

#include "../src/Face.h"

class ISurface {
  virtual double getDis(int index) = 0;

  virtual double getDis(int index, int cell) = 0;

  virtual std::optional<double> getFixed(int index) = 0;

  virtual const std::vector<Face> &getFaces() = 0;
};

#endif // ISURFACE_H
