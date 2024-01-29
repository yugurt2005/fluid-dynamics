#include "GridHelper.h"

Grid buildRectangularGrid(int n, int m, double s)
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

      centers[index] = {j + s / 2, i + s / 2};

      faces.push_back({index, getIndex(i, j + 1), Vector2d(j + 1, i) * s, Vector2d(j + 1, i + 1) * s});
      faces.push_back({getIndex(i + 1, j), index, Vector2d(j, i + 1) * s, Vector2d(j + 1, i + 1) * s});

      if (i == 0)
      {
        faces.push_back({index, -1, Vector2d(j, i) * s, Vector2d(j + 1, i) * s});
      }
      if (j == 0)
      {
        faces.push_back({-1, index, Vector2d(j, i) * s, Vector2d(j, i + 1) * s});
      }
    }
  }

  return Grid(centers, faces);
}
