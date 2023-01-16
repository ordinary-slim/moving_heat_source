#include <iostream>
#include <Eigen/Core>
#include <cmath>
#include <vector>
#include "includes/mesh.h"
#include "includes/element.h"
#include "includes/heatSource.h"

//1D heat sources
double gaussianPowerDensity(double x, double t, double x0, double power, double efficiency, double radius) {
  double pd = 2*(power*efficiency) / M_PI / pow(radius, 2) * exp( - 2*pow(x - x0, 2)/pow(radius, 2));
  return pd;
}

double forcedSolutionSource91(double x, double t, double x0, double cte, double gamma, double beta) {
  double alpha = 1;
  double pd = cte * exp( - gamma * t ) * exp( - beta * pow( x - x0, 2 ) );
  pd *= ( -pow(2 * beta * (x - x0), 2) + 2 * beta - gamma * alpha);
  return pd;
}

void HeatSource::computePulse( Eigen::VectorXd &pulse, Mesh &m, double t, double dt ) {
  double x_gp, r_i;
  updatePosition( dt );

  pulse.setZero();

  Element e;
  // assemble
  for (int ielem = 0; ielem < m.nels; ++ielem) {
    e = m.getElement( ielem );
    for (int inode = 0; inode < e.nnodes; ++inode) {
      r_i = 0;
      for (int igp = 0; igp < e.nnodes; ++igp) {
        x_gp = e.gpos[ igp ];
        r_i += e.gpweight[igp] * e.baseFunGpVals[inode][igp] * e.vol * powerDensity(x_gp, time, currentPosition[0], power, efficiency, radius);
      }
      pulse[e.con[inode]] += r_i;
    }
  }
}
