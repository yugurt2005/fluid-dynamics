#ifndef DEBUG_H
#define DEBUG_H

#include <cassert>
#include <iostream>
#include <iomanip>

#include "Eigen/Core"
#include "Eigen/Sparse"

using Eigen::VectorXd;

typedef Eigen::SparseMatrix<double> SpMat;

using std::cout;
using std::endl;

class Debug {
public:
  static void print(double input);

  static void debugSpMat(const SpMat &input);

  static void debug2d(const VectorXd &input, int n, int m, std::string name = "2d");

  static double debugResiduals(const VectorXd &old, const VectorXd &phi, std::string name);
};

#endif // DEBUG_H

