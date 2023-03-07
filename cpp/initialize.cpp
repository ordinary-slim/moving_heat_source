#include "problem.h"
#include <map>
#include <string>
#include <algorithm>
#include "../external/pybind11/include/pybind11/pybind11.h"
#include "../external/pybind11/include/pybind11/eigen.h"
#include "../external/pybind11/include/pybind11/stl.h"

namespace py = pybind11;

void Problem::initialize(py::dict &input) {
  // tstepping
  dt = py::cast<double>(input["dt"]);
  // MESH
  mesh.initializeMesh( input );

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
  mhs.currentPosition    = mhs.initialPosition;
  // set type of source term
  switch (int(py::cast<int>( input["sourceTerm"] ))) {
    case 91:
      { mhs.powerDensity = &forcedSolutionSource91;
        break; }
    default:
      {
        if (mesh.dim == 1 ) {
          mhs.powerDensity = &gaussianPowerDensity1D;
        } else if (mesh.dim == 2 ) {
          mhs.powerDensity = &gaussianPowerDensity2D;
        } else {
          printf("Dim > 2 not ready yet\n");
          exit(EXIT_FAILURE);
        }
        break; }
  }

  double environmentTemperature = py::cast<double>(input["environmentTemperature"]);
  // initialize unknown
  unknown.mesh = &mesh;
  unknown.values = Eigen::VectorXd::Constant( mesh.nnodes, environmentTemperature );

  // dirichlet BC
  vector<int> freeNodes(mesh.nnodes);


  // material. dictionnary is not efficient + involved in assembly
  material["rho"] = py::cast<double>(input["rho"]);
  material["k"] = py::cast<double>(input["conductivity"]);
  material["cp"] = py::cast<double>(input["specific_heat"]);

  // check for advection term
  if (input.contains("advectionSpeedX")) {
    advectionSpeed[0] = py::cast<double>(input["advectionSpeedX"]);
    advectionSpeed[1] = py::cast<double>(input["advectionSpeedY"]);
    advectionSpeed[2] = py::cast<double>(input["advectionSpeedZ"]);
    if (advectionSpeed.norm() > 1e-10) isAdvection = true;
  }

  // check for domain motion
  if (input.contains("speedFRF_X")) {
    mesh.speedFRF[0] = py::cast<double>(input["speedFRF_X"]);
    mesh.speedFRF[1] = py::cast<double>(input["speedFRF_Y"]);
    mesh.speedFRF[2] = py::cast<double>(input["speedFRF_Z"]);
  }

  // check for time dependency
  if (input.contains("steadyState")) {
    isSteady = py::cast<bool>(input["steadyState"]);
  }


  // timeIntegrator
  timeIntegrator.setRequiredSteps( py::cast<int>(input["timeIntegration"] ));

  // allocate storage for previous solutions
  unknown.prevValues = Eigen::MatrixXd::Zero( mesh.nnodes, timeIntegrator.nstepsRequired );
  unknown.prevValues.col( 0 ) << unknown.values;
  ++timeIntegrator.nstepsStored;
}
