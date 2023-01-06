#ifndef HEATSOURCE
#include <vector>
#include <Eigen/Core>
#include "mesh.h"

using namespace std;

class HeatSource {
  public:
      Eigen::Vector3d initialPosition;
      Eigen::Vector3d currentPosition;
      Eigen::Vector3d speed;
      double power;
      double efficiency = 1.0;
      double radius = 2.0;
      double time = 0.0;

      void updatePosition( double dt ) {
        currentPosition += speed * dt;
      }

      void computePulse( Eigen::VectorXd &pulse, Mesh &m, double t, double dt );
      double (*powerDensity)(double x, double t, double x0, double power, double efficiency, double radius);
    };

double gaussianPowerDensity(double x, double t, double x0, double power, double efficiency, double radius);
double forcedSolutionSource91(double x, double t, double x0, double power, double efficiency, double radius);
#define HEATSOURCE
#endif
