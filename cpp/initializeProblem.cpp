#include "Problem.h"
#include <map>
#include <string>
#include <algorithm>
#include "../external/pybind11/include/pybind11/pybind11.h"
#include "../external/pybind11/include/pybind11/eigen.h"
#include "../external/pybind11/include/pybind11/stl.h"

namespace py = pybind11;

Problem::Problem(mesh::Mesh &mesh, py::dict &input) :
  forcedDofs( mesh::MeshTag<int>( &mesh, 0, 0)),
  dirichletNodes( mesh::MeshTag<int>( &mesh ) ),
  dirichletValues( mesh::MeshTag<double>( &mesh ) ),
  neumannFacets( mesh::MeshTag<int>( &mesh, mesh.dim-1 ) ),
  neumannFluxes( mesh::MeshTag<vector<double>>( &mesh, mesh.dim-1 ) ),
  convectionFacets( mesh::MeshTag<int>( &mesh, mesh.dim-1 ) ),
  elsOwnedByOther( mesh::MeshTag<int>( &mesh, mesh.dim, 0 ) ),
  gammaNodes( mesh::MeshTag<int>( &mesh, 0 ) ),
  gammaFacets( mesh::MeshTag<int>( &mesh, mesh.dim-1 ) ),
  domain( mesh::ActiveMesh( &mesh ) ),
  unknown( fem::Function( &domain ) )
{

  // MATERIAL
  // TODO: Better DS!
  density = py::cast<double>(input["rho"]);
  conductivity = py::cast<double>(input["conductivity"]);
  specificHeat = py::cast<double>(input["specific_heat"]);
  if (input.contains("convectionCoeff")) {
    convectionCoeff = py::cast<double>(input["convectionCoeff"]);
    isConvection = true;
  }

  // HEAT SOURCE
  // set type of source term
  switch (int(py::cast<int>( input["sourceTerm"] ))) {
    case 86:
      { mhs = new cteHeat();
        break; }
    default:
      {
        if (domain.mesh->dim == 1 ) {
          mhs = new gaussianPowerDensity1D();
        } else if (domain.mesh->dim == 2 ) {
          mhs = new gaussianPowerDensity2D();
        } else {
          mhs = new gaussianPowerDensity3D();
        }
        break; }
  }
  mhs->radius = py::cast<double>(input["radius"]);
  mhs->power = py::cast<double>(input["power"]);
  mhs->pulse.resize( domain.mesh->nnodes );

  if (input.contains("efficiency")) mhs->efficiency = py::cast<double>(input["efficiency"]);

  mhs->speed[0] = py::cast<double>(input["HeatSourceSpeedX"]);
  mhs->speed[1] = py::cast<double>(input["HeatSourceSpeedY"]);
  mhs->speed[2] = py::cast<double>(input["HeatSourceSpeedZ"]);
  mhs->initialPosition[0] = py::cast<double>(input["initialPositionX"]);
  mhs->initialPosition[1] = py::cast<double>(input["initialPositionY"]);
  mhs->initialPosition[2] = py::cast<double>(input["initialPositionZ"]);
  mhs->currentPosition    = mhs->initialPosition;

  // TIME DEPENDENCY
  if (input.contains("steadyState")) {
    isSteady = py::cast<bool>(input["steadyState"]);
  }

  // INITIALIZE UNKNOWN
  Tenv = py::cast<double>(input["environmentTemperature"]);
  Tdeposition = Tenv;
  unknown.values = Eigen::VectorXd::Constant( domain.mesh->nnodes, Tenv );

  // TSTEPPING
  if ( not(isSteady)) {
    dt = py::cast<double>(input["dt"]);
    // TIME INTEGRATOR
    timeIntegrator.setRequiredSteps( py::cast<int>(input["timeIntegration"] ));
    // update time integrator
    previousValues.push_front(  unknown );
    ++timeIntegrator.nstepsStored;
  }

  // DIRICHLET BC
  if (input.contains("dirichletNodes")) {
    setDirichlet( py::cast<vector<int>>(input["dirichletNodes"]),
                  py::cast<vector<double>>(input["dirichletValues"]) );
  }

  // NEUMANN BC
  // TODO: Think about how to eat this

  // CONVECTION BC
  if (isConvection) {
    convectionFacets = domain.boundaryFacets;
  }

  // ADVECTION
  if (input.contains("advectionSpeedX")) {
    advectionSpeed[0] = py::cast<double>(input["advectionSpeedX"]);
    advectionSpeed[1] = py::cast<double>(input["advectionSpeedY"]);
    advectionSpeed[2] = py::cast<double>(input["advectionSpeedZ"]);
    if (advectionSpeed.norm() > 1e-10) isAdvection = true;
  }

  // DOMAIN MOTION
  // TODO: Move this down to mesh level!
  if (input.contains("speedFRF_X")) {
    domain.mesh->speedFRF[0] = py::cast<double>(input["speedFRF_X"]);
    domain.mesh->speedFRF[1] = py::cast<double>(input["speedFRF_Y"]);
    domain.mesh->speedFRF[2] = py::cast<double>(input["speedFRF_Z"]);
  }

  // ASSS STABILIZATION
  if (input.contains("isStabilized")) {
    isStabilized = py::cast<bool>(input["isStabilized"]);
  }
}
