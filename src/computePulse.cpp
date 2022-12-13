#include <iostream>
#include <Eigen/Core>
#include <cmath>
#include <vector>
#include "mesh.h"
#include "line.h"
#include "heatSource.h"

double HeatSource::powerDensity(double x) {
  double x0 = currentPosition[0];
  double pd = 2*(power*efficiency) / M_PI / pow(radius, 2) * exp( - 2*pow(x - x0, 2)/pow(radius, 2));
  return pd;
}

double HeatSource::ctePowerDensity(double t) {
  // For debugging purposes
  return (power*efficiency)*exp(-t);
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
        r_i += l.gpweight[igp] * l.baseFunGpVals[inode][igp] * l.vol * ctePowerDensity(t);
      }
      pulse[l.con[inode]] += r_i;
    }
  }
}
