#include "State.h"

State &State::relax(const State &old, double factor) {
  this->p = (1 - factor) * this->p + factor * old.p;
  this->k = (1 - factor) * this->k + factor * old.k;
  this->o = (1 - factor) * this->o + factor * old.o;

  // Return the current state
  return *this;
}