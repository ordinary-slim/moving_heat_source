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

vector<double> power_density( double x0, vector<double> x, double radius, double P ) {
  vector<double> pd( x.size() );
  for ( int idx = 0; idx < x.size(); ++idx ) {
    pd[idx] = 2*P / M_PI / pow(radius, 2) * exp( - 2*pow(x[idx] - x0, 2)/pow(radius, 2));
  }
  return pd;
}

int main() {
  double x0=0.0, radius=2, P=100, v = 10, L = 10, dt=0.1;
  int nels = 100, nSteps=100;
  Mesh m;
  m.initialize1DMesh(0, L, nels);

  for (int tstep=0; tstep<nSteps; ++tstep) {
    // update curr position
    x0 += tstep*dt*v;
    
    // get power density
    vector<double> pd = power_density(x0, m.pos, radius, P);

    plt::figure_size( 1200, 900 );
    plt::plot(m.pos, pd);

    plt::xlim(0.0, L);
    plt::title("Sample figure");
    //plt::legend();
    plt::show();
  }
}
