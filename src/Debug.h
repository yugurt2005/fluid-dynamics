#ifndef DEBUG_H
#define DEBUG_H

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
};

#endif // DEBUG_H

