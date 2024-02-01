#ifndef INLET_H
#define INLET_H

#include <optional>

using std::optional;

struct Inlet {
  optional<double> u;
  optional<double> v;
  optional<double> p;
};

#endif // INLET_H
