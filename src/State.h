#ifndef STATE_H
#define STATE_H

#include "Eigen/Core"

using Eigen::MatrixXd;

struct State {
  MatrixXd u;
  MatrixXd v;
  MatrixXd p;
  MatrixXd k;
  MatrixXd o;

  State &relax(const State &other, double factor);
};

#endif // STATE_H
