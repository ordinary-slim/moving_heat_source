#ifndef HEATSOURCE
#define HEATSOURCE
#include <vector>
#include <Eigen/Core>
#include <cmath>

namespace heat {
enum HeatSourceType {
  none,
  gaussian1d,
  gaussian2d,
  gaussian3d,
  constant,
  lumped,
};

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
      HeatSourceType type = none;

      void updatePosition( double dt ) {
        currentPosition += speed * dt;
      }
      void setSpeed( Eigen::Vector3d inputSpeed ) {
        speed = inputSpeed;
      }
      void setPower( double otherPower ) {
        power = otherPower;
      }

      //double (*powerDensity)(Eigen::Vector3d x, double t);

      virtual double operator()(Eigen::Vector3d x, double t) const {
        return 0.0;
      }
};

// Some heat sources
class gaussianPowerDensity1D : public HeatSource {
  public:
    gaussianPowerDensity1D() {
      type = gaussian1d;
    }
    double operator()(Eigen::Vector3d x, double t) const {
      double pd = 2*(power*efficiency) / M_PI / pow(radius, 2) * exp( - 2*pow(x[0] - currentPosition[0], 2)/pow(radius, 2));
      return pd;
    }
};

class gaussianPowerDensity2D : public HeatSource {
  public:
    gaussianPowerDensity2D() {
      type = gaussian2d;
    }
    double operator()(Eigen::Vector3d x, double t) const {
      double pd = 2*(power*efficiency) / M_PI / pow(radius, 2) * exp( - 2*(pow(x[0] - currentPosition[0], 2) + pow(x[1] - currentPosition[1], 2)) /pow(radius, 2));
      return pd;
    }
};

class gaussianPowerDensity3D : public HeatSource {
  public:
    gaussianPowerDensity3D() {
      type = gaussian3d;
    }
    double operator()(Eigen::Vector3d x, double t) const {
      //TODO: check formula
      double dSquared = (x - currentPosition).squaredNorm();
      double pd = 6*sqrt(3)*(power*efficiency) / pow(M_PI, 1.5) / pow(radius, 3) * exp( -3*dSquared/pow(radius, 2));
      return pd;
    }
};

class cteHeat : public HeatSource {
  public:
    cteHeat() {
      type = constant;
    }
    double operator()(Eigen::Vector3d x, double t) const {

      double pd = 0;
      if ( (x - currentPosition).norm() <= radius ) {
        //pd = efficiency * power / 2 / radius;
        pd = power;
      }
      return pd;
    }
};
}
#endif
