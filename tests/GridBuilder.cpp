#include "GridBuilder.h"

vector<Face> GridBuilder::buildRectangularGrid(int n, int m, double w, double h)
{
  vector<Face> faces;

  auto getIndex = [n, m](int i, int j)
  {
    if (i < 0 || i >= n || j < 0 || j >= m)
      return -1;
    return i * m + j;
  };

  for (int i = 0; i <= n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      int l = getIndex(i - 0, j);
      int r = getIndex(i - 1, j);
      faces.push_back({l, r, {(j + 1) * w, i * h}, {j * w, i * h}});
    }
  }

  for (int j = 0; j <= m; j++)
  {
    for (int i = 0; i < n; i++)
    {
      int l = getIndex(i, j - 1);
      int r = getIndex(i, j - 0);
      faces.push_back({l, r, {j * w, (i + 1) * h}, {j * w, i * h}});
    }
  }

  return faces;
}

vector<Face> GridBuilder::buildTriangularGrid(int n, int m, double w, double h)
{
  vector<Face> faces;

  auto getIndex = [n, m](int i, int j, bool side)
  {
    if (i < 0 || i >= n || j < 0 || j >= m)
      return -1;
    return (i * m + j) * 2 + side;
  };

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      int l = getIndex(i, j, 0);
      int r = getIndex(i, j, 1);
      faces.push_back({l, r, {j * w, (i + 1) * h}, {(j + 1) * w, i * h}});
    }
  }

  for (int i = 0; i <= n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      int l = getIndex(i - 0, j, 0);
      int r = getIndex(i - 1, j, 1);
      faces.push_back({l, r, {(j + 1) * w, i * h}, {j * w, i * h}});
    }
  }

  for (int j = 0; j <= m; j++)
  {
    for (int i = 0; i < n; i++)
    {
      int l = getIndex(i, j - 1, 1);
      int r = getIndex(i, j - 0, 0);
      faces.push_back({l, r, {j * w, (i + 1) * h}, {j * w, i * h}});
    }
  }

  return faces;
}