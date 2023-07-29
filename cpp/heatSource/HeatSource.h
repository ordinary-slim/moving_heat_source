#ifndef HEATSOURCE
#define HEATSOURCE
#include <vector>
#include <Eigen/Core>
#include <cmath>
#include <memory>
#include "Path.h"

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
      std::unique_ptr<Path> path = NULL;
      Eigen::Vector3d initialPosition;
      Eigen::Vector3d currentPosition;
      Eigen::Vector3d speed;
      Eigen::VectorXd pulse; // source term
      double power;
      double efficiency = 1.0;
      double radius = 2.0;
      double time = 0.0;// Is this used ?
      HeatSourceType type = none;

      void updatePosition( double dt ) { currentPosition += speed * dt; }
      void setSpeed( Eigen::Vector3d speed ) { this->speed = speed; }
      void setPower( double power ) { this->power = power; }
      void setPath( std::vector<Eigen::Vector3d> &coordinates,
            std::vector<double> &speeds,
            std::vector<double> &powers,
            std::vector<int> &arePrinting ) {

          this->path = std::make_unique<heat::Path>( coordinates, speeds, powers, arePrinting );
      }
      void preIterate(double newTime) {
        double dt = newTime - this->time;
        this->time = newTime;
        if (path != NULL) {
          path->currentTrack = path->interpolateTrack( this->time );
          if (path->currentTrack != NULL) {
            this->speed = path->currentTrack->getSpeed();
            this->power = path->currentTrack->power;
          } else {
            throw std::invalid_argument("Time is out of bounds.");
          }
        }
        updatePosition( dt );
      }

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
