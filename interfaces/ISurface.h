#ifndef ISURFACE_H
#define ISURFACE_H

#include <optional>

#include "../src/Face.h"

class ISurface {
public:
  virtual double getDis(int index) const = 0;

  virtual double getDis(int index, int cell) const = 0;

  virtual std::optional<double> getFixed(int index) const = 0;

  virtual const std::vector<Face> &getFaces() const = 0;
};

#endif // ISURFACE_H
