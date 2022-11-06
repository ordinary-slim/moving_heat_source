#include <iostream>
#include <Eigen/Core>
#include <cmath>
#include <matplotlib-cpp/matplotlibcpp.h>
#include <vector>
#include "mesh.h"
#include "line.h"

namespace plt = matplotlibcpp;

vector<double> compute_pulse(double x, double radius, double P, Mesh m) {
  vector<double> r( m.nnodes );

  Line l;
  for (int ielem = 0; ielem < m.nels; ++ielem) {
    for (int igp = 0; igp < l.nnodes; ++igp) {
    }
  }

  return r;
}

double power_density( double x0, double x, double radius, double P ) {
  double pd = 2*P / M_PI / pow(radius, 2) * exp( - 2*pow(x - x0, 2)/pow(radius, 2));
  return pd;
}

int main() {
  double x = 4, x0=5, radius=2, P=100;
  cout << "pd:" << power_density( x0, x, radius, P ) << endl;
}
