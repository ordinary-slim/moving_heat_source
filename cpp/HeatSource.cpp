#include <iostream>
#include <Eigen/Core>
#include <cmath>
#include "HeatSource.h"

double gaussianPowerDensity1D(Eigen::Vector3d x, double t, Eigen::Vector3d x0, double power, double efficiency, double radius) {
  double pd = 2*(power*efficiency) / M_PI / pow(radius, 2) * exp( - 2*pow(x[0] - x0[0], 2)/pow(radius, 2));
  return pd;
}
double gaussianPowerDensity2D(Eigen::Vector3d x, double t, Eigen::Vector3d x0, double power, double efficiency, double radius) {
  double pd = 2*(power*efficiency) / M_PI / pow(radius, 2) * exp( - 2*(pow(x[0] - x0[0], 2) + pow(x[1] - x0[1], 2)) /pow(radius, 2));
  return pd;
}

double gaussianPowerDensity3D(Eigen::Vector3d x, double t, Eigen::Vector3d x0, double power, double efficiency, double radius) {
  //TODO: check formula
  double dSquared = (x - x0).squaredNorm();
  double pd = 6*sqrt(3)*(power*efficiency) / pow(M_PI, 1.5) / pow(radius, 3) * exp( -3*dSquared/pow(radius, 2));
  return pd;
}

double cteHeat(Eigen::Vector3d x, double t, Eigen::Vector3d x0,
    double power, double efficiency, double radius) {

  double pd = 0;
  if ( (x - x0).norm() <= radius ) {
    //pd = efficiency * power / 2 / radius;
    pd = power;
  }
  return pd;
}

double forcedSolutionSource91(Eigen::Vector3d x, double t, Eigen::Vector3d x0, double cte, double gamma, double beta) {
  double alpha = 1;
  double pd = cte * exp( - gamma * t ) * exp( - beta * pow( x[0] - x0[0], 2 ) );
  pd *= ( -pow(2 * beta * (x[0] - x0[0]), 2) + 2 * beta - gamma * alpha);
  return pd;
}
