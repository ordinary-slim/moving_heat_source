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

      void updatePosition( double t ) {
        currentPosition = initialPosition + speed * t;
      }

      void computePulse( Eigen::VectorXd &pulse, double t, Mesh &m );
      double powerDensity(double x);
      double ctePowerDensity(double t);
    };
#define HEATSOURCE
#endif
