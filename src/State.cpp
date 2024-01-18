#include "State.h"

State &State::relax(const State &other, double factor) {
  // Modify the current state based on the other state and the factor
  // For example, if you want to average the current and other states:
  this->u = (1 - factor) * this->u + factor * other.u;
  this->v = (1 - factor) * this->v + factor * other.v;
  this->p = (1 - factor) * this->p + factor * other.p;
  this->k = (1 - factor) * this->k + factor * other.k;
  this->o = (1 - factor) * this->o + factor * other.o;

  // Return the current state
  return *this;
}