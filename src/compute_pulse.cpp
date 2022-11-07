#include <iostream>
#include <Eigen/Core>
#include <cmath>
#include <vector>
#include "mesh.h"
#include "line.h"

#include <matplotlib-cpp/matplotlibcpp.h>
namespace plt = matplotlibcpp;


double power_density( double x0, double x, double radius, double P ) {
  double pd = 2*P / M_PI / pow(radius, 2) * exp( - 2*pow(x - x0, 2)/pow(radius, 2));
  return pd;
}

void compute_pulse(Eigen::VectorXd &r, double x0, double radius, double P, Mesh &m) {
  double x_gp, r_i;
  Line l;
  for (int ielem = 0; ielem < m.nels; ++ielem) {
    l = m.getElement( ielem );
    for (int inode = 0; inode < l.nnodes; ++inode) {
      r_i = 0;
      for (int igp = 0; igp < l.nnodes; ++igp) {
        x_gp = l.gpos[ igp ];
        r_i += l.gpweight[igp] * l.baseFunGpVals[inode][igp] * power_density( x0, x_gp, radius, P);
      }
      r[inode] += r_i;
    }
  }
}

/*
int main() {
double x0=0.0, radius=2, P=100, v = 10, L = 10, dt=0.05, maxPower = -1;
int nels = 100, nSteps=10;
Mesh m;
m.initialize1DMesh(0, L, nels);
vector<double> r( m.nnodes ); // source term

for (int tstep=0; tstep<nSteps; ++tstep) {
  // update curr position
  x0 += tstep*dt*v;
  // get power density
  compute_pulse(r, x0, radius, P, m);
  
  plt::clf();
  plt::plot( m.pos, r );
  plt::xlim( 0.0, L );
  double maxPower = max( maxPower, *max_element( r.begin(), r.end() ) );
  plt::ylim( 0.0, 1.1*maxPower );

  plt::pause( 0.04 );
}
}
*/
