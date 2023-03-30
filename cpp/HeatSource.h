#ifndef HEATSOURCE
#include <vector>
#include <Eigen/Core>

class HeatSource {
  public:
      Eigen::Vector3d initialPosition;
      Eigen::Vector3d currentPosition;
      Eigen::Vector3d speed;
      Eigen::VectorXd pulse; // source term
      double power;
      double efficiency = 1.0;
      double radius = 2.0;
      double time = 0.0;

      void updatePosition( double dt ) {
        currentPosition += speed * dt;
      }
      void setSpeed( Eigen::Vector3d inputSpeed ) {
        speed = inputSpeed;
      }
      void setPower( double otherPower ) {
        power = otherPower;
      }

      double (*powerDensity)(Eigen::Vector3d x, double t, Eigen::Vector3d x0, double power, double efficiency, double radius);
    };

// Some heat sources
double gaussianPowerDensity1D(Eigen::Vector3d x, double t,
    Eigen::Vector3d x0, double power, double efficiency, double radius);

double gaussianPowerDensity2D(Eigen::Vector3d x, double t,
    Eigen::Vector3d x0, double power, double efficiency, double radius);

double testCteHeat1D(Eigen::Vector3d x, double t, Eigen::Vector3d x0,
    double power, double efficiency, double radius);

double forcedSolutionSource91(Eigen::Vector3d x, double t,
    Eigen::Vector3d x0, double power, double efficiency, double radius);
#define HEATSOURCE
#endif
