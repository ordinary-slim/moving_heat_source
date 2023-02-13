#include <iostream>
#include <Eigen/Core>
#include <cmath>
#include <vector>
#include "domain/mesh.h"
#include "domain/element.h"
#include "heatSource.h"

//1D heat sources
double gaussianPowerDensity1D(Eigen::Vector3d x, double t, Eigen::Vector3d x0, double power, double efficiency, double radius) {
  double pd = 2*(power*efficiency) / M_PI / pow(radius, 2) * exp( - 2*pow(x[0] - x0[0], 2)/pow(radius, 2));
  return pd;
}
double gaussianPowerDensity2D(Eigen::Vector3d x, double t, Eigen::Vector3d x0, double power, double efficiency, double radius) {
  //TODO: check formula
  double pd = 2*(power*efficiency) / M_PI / pow(radius, 2) * exp( - 2*(pow(x[0] - x0[0], 2) + pow(x[1] - x0[1], 2)) /pow(radius, 2));
  return pd;
}

double forcedSolutionSource91(Eigen::Vector3d x, double t, Eigen::Vector3d x0, double cte, double gamma, double beta) {
  double alpha = 1;
  double pd = cte * exp( - gamma * t ) * exp( - beta * pow( x[0] - x0[0], 2 ) );
  pd *= ( -pow(2 * beta * (x[0] - x0[0]), 2) + 2 * beta - gamma * alpha);
  return pd;
}

void HeatSource::computePulse( Eigen::VectorXd &pulse, Mesh &m, double t, double dt ) {
  double r_i;
  Eigen::Vector3d x_gp;

  updatePosition( dt );

  pulse.setZero();

  Element e;
  // assemble
  for (int ielem = 0; ielem < m.nels; ++ielem) {
    if (m.activeElements[ielem]==0){
      continue;
    }
    e = m.getElement( ielem );
    for (int inode = 0; inode < e.nnodes; ++inode) {
      r_i = 0;
      for (int igp = 0; igp < e.nnodes; ++igp) {
        x_gp = e.gpos.row( igp );
        r_i += e.gpweight[igp] * e.BaseGpVals[inode][igp] * e.vol * powerDensity(x_gp, time, currentPosition, power, efficiency, radius);
      }
      pulse[e.con[inode]] += r_i;
    }
  }
}
