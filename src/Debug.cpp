#include "Debug.h"

void Debug::print(double input)
{
  cout << std::setw(3) << input;
}

void Debug::debugSpMat(const SpMat &input)
{
  int n = input.rows();
  int m = input.cols();

  cout << "SpMat: " << n << "x" << m << endl;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      print(input.coeff(i, j));
    }
    cout << endl;
  }
}