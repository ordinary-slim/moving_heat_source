#include <iostream>
#include <Eigen/Core>
#include <cmath>
#include <vector>
#include "includes/mesh.h"
#include "includes/line.h"
#include "includes/heatSource.h"

double gaussianPowerDensity(double x, double x0, double power, double efficiency, double radius) {
  double pd = 2*(power*efficiency) / M_PI / pow(radius, 2) * exp( - 2*pow(x - x0, 2)/pow(radius, 2));
  return pd;
}

double gaussianPowerDensityMRF(double x, double x0, double power, double efficiency, double radius) {
  double pd = 2*(power*efficiency) / M_PI / pow(radius, 2) * exp( - 2*pow(x, 2)/pow(radius, 2));
  return pd;
}

void HeatSource::computePulse( Eigen::VectorXd &pulse, double t, Mesh &m ) {
  double x_gp, r_i;
  updatePosition( t );

  pulse.setZero();

  Line l;
  // assemble
  for (int ielem = 0; ielem < m.nels; ++ielem) {
    l = m.getElement( ielem );
    for (int inode = 0; inode < l.nnodes; ++inode) {
      r_i = 0;
      for (int igp = 0; igp < l.nnodes; ++igp) {
        x_gp = l.gpos[ igp ];
        r_i += l.gpweight[igp] * l.baseFunGpVals[inode][igp] * l.vol * powerDensity(x_gp, currentPosition[0], power, efficiency, radius);
      }
      pulse[l.con[inode]] += r_i;
    }
  }
}
