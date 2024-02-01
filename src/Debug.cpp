#include "Debug.h"

void Debug::print(double input)
{
  cout << std::setw(4) << std::setprecision(2) << std::fixed << input;
}

void Debug::debugSpMat(const SpMat &input)
{
  int n = input.rows();
  int m = input.cols();

  cout << "SpMat: " << n << "x" << m << endl;
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      print(input.coeff(i, j));
      cout << " ";
    }
    cout << endl;
  }
}

void Debug::debug2d(const VectorXd &input, int n, int m, std::string name)
{
  assert(input.size() == n * m);

  cout << name << ": " << n << "x" << m << endl;
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      print(input(i));
      cout << " ";
    }
    cout << endl;
  }
  cout << endl;
}

double Debug::debugResiduals(const VectorXd &phi, const VectorXd &old, std::string name)
{
  double residual = (phi - old).norm();
  cout << name << " Residual: " << residual << endl << endl;
  return residual;
}