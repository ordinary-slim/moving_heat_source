#include "Problem.h"
#include <map>
#include <string>
#include <algorithm>
#include "../external/pybind11/include/pybind11/pybind11.h"
#include "../external/pybind11/include/pybind11/eigen.h"
#include "../external/pybind11/include/pybind11/stl.h"

namespace py = pybind11;

void Problem::initialize(py::dict &input) {
  // MESH
  mesh.initializeMesh( input );

  // MATERIAL
  // TODO: Better DS!
  material["rho"] = py::cast<double>(input["rho"]);
  material["k"] = py::cast<double>(input["conductivity"]);
  material["cp"] = py::cast<double>(input["specific_heat"]);
  if (input.contains("convectionCoeff")) {
    material["h"] = py::cast<double>(input["convectionCoeff"]);
    isConvection = true;
  }

  // HEAT SOURCE
  mhs.radius = py::cast<double>(input["radius"]);
  mhs.power = py::cast<double>(input["power"]);
  mhs.pulse.resize( mesh.nnodes );

  if (input.contains("efficiency")) mhs.efficiency = py::cast<double>(input["efficiency"]);

  mhs.speed[0] = py::cast<double>(input["HeatSourceSpeedX"]);
  mhs.speed[1] = py::cast<double>(input["HeatSourceSpeedY"]);
  mhs.speed[2] = py::cast<double>(input["HeatSourceSpeedZ"]);
  mhs.initialPosition[0] = py::cast<double>(input["initialPositionX"]);
  mhs.initialPosition[1] = py::cast<double>(input["initialPositionY"]);
  mhs.initialPosition[2] = py::cast<double>(input["initialPositionZ"]);
  mhs.currentPosition    = mhs.initialPosition;
  // set type of source term
  switch (int(py::cast<int>( input["sourceTerm"] ))) {
    case 91:
      { mhs.powerDensity = &forcedSolutionSource91;
        break; }
    case 86:
      { mhs.powerDensity = &testCteHeat1D;
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

  // TIME DEPENDENCY
  if (input.contains("steadyState")) {
    isSteady = py::cast<bool>(input["steadyState"]);
  } else {
    // TSTEPPING
    dt = py::cast<double>(input["dt"]);
  }

  // TIME INTEGRATOR
  timeIntegrator.setRequiredSteps( py::cast<int>(input["timeIntegration"] ));
  // INITIALIZE UNKNOWN
  Tenv = py::cast<double>(input["environmentTemperature"]);
  unknown = FEMFunction( mesh, timeIntegrator.nstepsRequired );
  unknown.values = Eigen::VectorXd::Constant( mesh.nnodes, Tenv );
  // update time integrator
  unknown.prevValues.col( 0 ) << unknown.values;
  ++timeIntegrator.nstepsStored;



  // DIRICHLET BC
  if (input.contains("dirichletNodes")) {
    unknown.dirichletNodes = py::cast<vector<int>>(input["dirichletNodes"]);
    unknown.dirichletValues = py::cast<vector<double>>(input["dirichletValues"]);
  }

  // NEUMANN BC
  // TODO: Think about how to eat this

  // CONVECTION BC
  if (input.contains("convectionCoeff")) {
    convectionFacets = mesh.boundaryFacets;
  }

  // ADVECTION
  if (input.contains("advectionSpeedX")) {
    advectionSpeed[0] = py::cast<double>(input["advectionSpeedX"]);
    advectionSpeed[1] = py::cast<double>(input["advectionSpeedY"]);
    advectionSpeed[2] = py::cast<double>(input["advectionSpeedZ"]);
    if (advectionSpeed.norm() > 1e-10) isAdvection = true;
  }

  // DOMAIN MOTION
  if (input.contains("speedFRF_X")) {
    mesh.speedFRF[0] = py::cast<double>(input["speedFRF_X"]);
    mesh.speedFRF[1] = py::cast<double>(input["speedFRF_Y"]);
    mesh.speedFRF[2] = py::cast<double>(input["speedFRF_Z"]);
  }

  // ASSS STABILIZATION
  if (input.contains("isStabilized")) {
    isStabilized = py::cast<bool>(input["isStabilized"]);
  }

}
