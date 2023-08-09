#ifndef HEATSOURCE
#define HEATSOURCE
#include <vector>
#include <Eigen/Core>
#include <cmath>
#include <memory>
#include "Path.h"
#include "../../external/pybind11/include/pybind11/pybind11.h"
#include "../../external/pybind11/include/pybind11/numpy.h"
#include "../../external/pybind11/include/pybind11/eigen.h"
#include "../../external/pybind11/include/pybind11/stl.h"

class Problem;

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
  protected:
      const Problem *problem;//observer
  public:
      std::unique_ptr<Path> path = NULL;
      Eigen::Vector3d position;
      Eigen::Vector3d speed;
      Eigen::VectorXd pulse; // source term
      double power;
      double efficiency = 1.0;
      double radius = 2.0;
      HeatSourceType type = none;
      const Track *currentTrack = NULL;

      HeatSource( pybind11::dict &input, Problem *problem );

      void step( double dt ) { position += speed * dt; }
      void setSpeed( Eigen::Vector3d speed ) { this->speed = speed; }
      void setPower( double power ) { this->power = power; }
      void setPath( std::vector<Eigen::Vector3d> &coordinates,
            std::vector<double> &speeds,
            std::vector<double> &powers,
            std::vector<int> &arePrinting ) {

          path = std::make_unique<heat::Path>( coordinates, speeds, powers, arePrinting );
          currentTrack = &path->tracks[0];
      }

      void preIterate();
      virtual double operator()(Eigen::Vector3d x, double t) const {
        return 0.0;
      }
};

// Some heat sources
class gaussianPowerDensity1D : public HeatSource {
  public:
    gaussianPowerDensity1D( pybind11::dict &input, Problem *problem )
      : HeatSource(input, problem)
    {
      type = gaussian1d;
    }
    double operator()(Eigen::Vector3d x, double t) const {
      double pd = 2*(power*efficiency) / M_PI / pow(radius, 2) * exp( - 2*pow(x[0] - position[0], 2)/pow(radius, 2));
      return pd;
    }
};

class gaussianPowerDensity2D : public HeatSource {
  public:
    gaussianPowerDensity2D( pybind11::dict &input, Problem *problem )
      : HeatSource(input, problem)
    {
      type = gaussian2d;
    }
    double operator()(Eigen::Vector3d x, double t) const {
      double pd = 2*(power*efficiency) / M_PI / pow(radius, 2) * exp( - 2*(pow(x[0] - position[0], 2) + pow(x[1] - position[1], 2)) /pow(radius, 2));
      return pd;
    }
};

class gaussianPowerDensity3D : public HeatSource {
  public:
    gaussianPowerDensity3D( pybind11::dict &input, Problem *problem )
      : HeatSource(input, problem)
    {
      type = gaussian3d;
    }
    double operator()(Eigen::Vector3d x, double t) const {
      //TODO: check formula
      double dSquared = (x - position).squaredNorm();
      double pd = 6*sqrt(3)*(power*efficiency) / pow(M_PI, 1.5) / pow(radius, 3) * exp( -3*dSquared/pow(radius, 2));
      return pd;
    }
};

class cteHeat : public HeatSource {
  public:
    cteHeat( pybind11::dict &input, Problem *problem )
      : HeatSource(input, problem)
    {
      type = constant;
    }
    double operator()(Eigen::Vector3d x, double t) const {

      double pd = 0;
      if ( (x - position).norm() <= radius ) {
        //pd = efficiency * power / 2 / radius;
        pd = power;
      }
      return pd;
    }
};
}
#endif
