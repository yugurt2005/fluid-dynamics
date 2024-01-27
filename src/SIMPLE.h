#ifndef SIMPLE_H
#define SIMPLE_H

#include "Eigen/Dense"
#include "Eigen/IterativeLinearSolvers"
#include "Eigen/Sparse"

#include "State.h"
#include "Parameters.h"
#include "FVM.h"

using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::Vector2d;
using Eigen::VectorXd;

typedef SparseMatrix<double> SpMat;

class SIMPLE
{
};

#endif // SIMPLE_H