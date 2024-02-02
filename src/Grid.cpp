#include "Grid.h"

Grid::Grid(const vector<Face> &faces) {
  m = faces.size();
  n = 0;

  for (const Face &f : faces) {
    n = std::max(m, f.l);
    n = std::max(m, f.r);
  }

  volume = new double[n];
  center = new Vector[n];

  edges = new vector<Edge>[n]();
  for (int i = 0; i < m; i++) {
    const Face &f = faces[i];
    if (f.l != -1) {
      edges[f.l].push_back(Edge(f.r, i, f));
    }
    if (f.r != -1) {
      edges[f.r].push_back(Edge(f.l, i, f));
    }
  }

  for (int i = 0; i < n; i++) {
    volume[i] = calcVolume(edges[i]);
    center[i] = calcCenter(edges[i], volume[i]);
  }
}

Grid::~Grid() {
  delete[] volume;
  delete[] center;
  delete[] edges;
}

double Grid::calcVolume(const vector<Edge> &edges) {
  double volume = 0;
  for (const Edge &e : edges) {
    volume += e.a * e.c.dot(e.n) / 2;
  }
  return volume;
}

Vector Grid::calcCenter(const vector<Edge> &edges, double volume) {
  Vector t(0, 0);
  for (const Edge &e : edges) {
    t += e.c;
  }
  t /= edges.size();

  Vector c(0, 0);
  for (const Edge &e : edges) {
    Vector tp = e.p - t;
    Vector tq = e.q - t;
    double triangle = std::abs(tp.y() * tq.x() - tp.x() * tq.y()) / 2;

    c += (e.p + e.q + t) / 3 * triangle / volume;
  }

  return c;
}