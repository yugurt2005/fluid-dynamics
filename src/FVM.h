#ifndef CFD_FVM_H
#define CFD_FVM_H

#include <vector>

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/Sparse"

#include "../interfaces/IGrid.h"

using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::Vector2d;
using Eigen::VectorXd;
using std::vector;

typedef SparseMatrix<double> SpMat;

class FVM
{
private:
  int n;
  int m;

  IGrid &grid;

  SparseMatrix<double> gradientX;
  SparseMatrix<double> gradientY;

  inline vector<int> getAdjacents(int index) { return grid.getAdjacents(index); }

  inline VectorXd getAdjDx(int index) { return grid.getAdjDx(index); }

  inline VectorXd getAdjDy(int index) { return grid.getAdjDy(index); };

  inline vector<Face> &getFaces() { return grid.getFaces(); }

public:
  FVM(IGrid &grid);

  void init();

  static VectorXd calcLeastSquaresGradient(int n, const VectorXd &dx, const VectorXd &dy);

  void createGradientMatrix();

  inline const SpMat &getGradientX() { return gradientX; };

  inline const SpMat &getGradientY() { return gradientY; };

  double interpolateAt(int index, Face &face, double phi, double dx, double dy);

  VectorXd interpolate(VectorXd &phi, VectorXd &flux);

  double calcMassFluxAt(int index, Face &face, Vector2d velocity);

  VectorXd calcMassFlux(const VectorXd &u, const VectorXd &v);
};

#endif // CFD_FVM_H