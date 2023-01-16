#include "includes/problem.h"
#include <map>
#include <string>
#include "../external/pybind11/include/pybind11/pybind11.h"
#include "../external/pybind11/include/pybind11/stl.h"

namespace py = pybind11;

void Problem::initialize(py::dict &input) {
  // tstepping
  dt = py::cast<double>(input["dt"]);
  // discrete mesh
  if (input.contains("1D")) {
    double a = py::cast<double>(input["Left"]);
    double b = py::cast<double>(input["Right"]);
    int nels = py::cast<int>( input["nels"] );
    mesh.initialize1DMesh( a, b, nels );
  }

  // heat source
  mhs.radius = py::cast<double>(input["radius"]);
  mhs.power = py::cast<double>(input["power"]);
  if (input.contains("efficiency")) mhs.efficiency = py::cast<double>(input["efficiency"]);

  mhs.speed[0] = py::cast<double>(input["speedX"]);
  mhs.speed[1] = py::cast<double>(input["speedY"]);
  mhs.speed[2] = py::cast<double>(input["speedZ"]);
  mhs.initialPosition[0] = py::cast<double>(input["initialPositionX"]);
  mhs.initialPosition[1] = py::cast<double>(input["initialPositionY"]);
  mhs.initialPosition[2] = py::cast<double>(input["initialPositionZ"]);
  // set type of source term
  switch (int(py::cast<int>( input["sourceTerm"] ))) {
    case 91:
      { mhs.powerDensity = &forcedSolutionSource91;
        break; }
    default:
      { mhs.powerDensity = &gaussianPowerDensity;
        break; }
  }

  // initialize solution and increment
  double environmentTemperature = py::cast<double>(input["environmentTemperature"]);
  solution = Eigen::VectorXd::Constant( mesh.nnodes, environmentTemperature );

  // dirichlet BC
  vector<int> freeNodes(mesh.nnodes);


  // material. dictionnary is not efficient + involved in assembly
  material["rho"] = py::cast<double>(input["rho"]);
  material["k"] = py::cast<double>(input["conductivity"]);
  material["cp"] = py::cast<double>(input["specific_heat"]);

  // check for advection term
  if (input.contains("isAdvection")) {
    isAdvection = py::cast<bool>(input["isAdvection"]);
    if (isAdvection) {
      advectionSpeed[0] = py::cast<double>(input["advectionSpeedX"]);
      advectionSpeed[1] = py::cast<double>(input["advectionSpeedY"]);
      advectionSpeed[2] = py::cast<double>(input["advectionSpeedZ"]);
      cout << "advectionSpeed= " << advectionSpeed << endl;
    }
  }

  // check for time dependency
  if (input.contains("steadyState")) {
    isSteady = py::cast<bool>(input["steadyState"]);
  }


  // timeIntegrator
  timeIntegrator.setRequiredSteps( py::cast<int>(input["timeIntegration"] ));

  // allocate storage for previous solutions
  prevSolutions = Eigen::MatrixXd::Zero( mesh.nnodes, timeIntegrator.nstepsRequired );
  prevSolutions.col( 0 ) << solution;
  ++timeIntegrator.nstepsStored;
}
