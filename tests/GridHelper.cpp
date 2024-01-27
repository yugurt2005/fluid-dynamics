#include "GridHelper.h"

Grid buildRectangularGrid(int n, int m)
{
  auto getIndex = [n, m](int i, int j)
  {
    return i < 0 || i >= n || j < 0 || j >= m ? -1 : i * m + j;
  };

  vector<Vector2d> centers(n * m);
  vector<Face> faces;

  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      int index = getIndex(i, j);

      assert(index < n * m);
      assert(index >= 0);

      centers[index] = {j + 0.5, i + 0.5};

      faces.push_back({index, getIndex(i, j + 1), {j + 1, i}, {j + 1, i + 1}});
      faces.push_back({getIndex(i + 1, j), index, {j, i + 1}, {j + 1, i + 1}});

      if (i == 0)
      {
        faces.push_back({index, -1, {j, i}, {j + 1, i}});
      }
      if (j == 0)
      {
        faces.push_back({-1, index, {j, i}, {j, i + 1}});
      }
    }
  }

  return Grid(centers, faces);
}
