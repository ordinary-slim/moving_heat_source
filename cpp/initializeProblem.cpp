#include "Problem.h"
#include <map>
#include <string>
#include <algorithm>
#include "../external/pybind11/include/pybind11/pybind11.h"
#include "../external/pybind11/include/pybind11/numpy.h"
#include "../external/pybind11/include/pybind11/eigen.h"
#include "../external/pybind11/include/pybind11/stl.h"

namespace py = pybind11;

Eigen::VectorXd CreateEigenVector(py::array_t<double> n);

Problem::Problem(mesh::Mesh &mesh, py::dict &input) :
  domain( mesh::ActiveMesh( &mesh ) ),
  unknown( fem::Function( &domain ) ),
  forcedDofs( mesh::MeshTag<int>( &mesh, 0, 0)),
  dirichletNodes( mesh::MeshTag<int>( &mesh ) ),
  dirichletValues( mesh::MeshTag<double>( &mesh ) ),
  weakBcFacets( mesh::MeshTag<int>( &mesh, mesh.dim-1 ) ),
  neumannFluxes( mesh::MeshTag<vector<double>>( &mesh, mesh.dim-1, vector<double>() ) ),
  elsOwnedByOther( mesh::MeshTag<int>( &mesh, mesh.dim, 0 ) ),
  gammaNodes( mesh::MeshTag<int>( &mesh, 0 ) ),
  gammaFacets( mesh::MeshTag<int>( &mesh, mesh.dim-1 ) )
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
    case 11: { 
        double heatSouceWidth  = py::cast<double>(input["heatSourceWidth"]);
        double heatSouceHeight = py::cast<double>(input["heatSourceHeight"]);
        mhs = std::make_unique<heat::LumpedHeatSource>( &domain, heatSouceWidth, heatSouceHeight);
        break; }
    case 86: {
        mhs = std::make_unique<heat::cteHeat>();
        break; }
    default:
      {
        if (domain.mesh->dim == 1 ) {
          mhs = std::make_unique<heat::gaussianPowerDensity1D>();
        } else if (domain.mesh->dim == 2 ) {
          mhs = std::make_unique<heat::gaussianPowerDensity2D>();
        } else {
          mhs = std::make_unique<heat::gaussianPowerDensity3D>();
        }
        break; }
  }
  mhs->radius = py::cast<double>(input["radius"]);
  mhs->power = py::cast<double>(input["power"]);
  mhs->pulse.resize( domain.mesh->nnodes );

  if (input.contains("efficiency")) mhs->efficiency = py::cast<double>(input["efficiency"]);

  mhs->speed = CreateEigenVector(py::array_t<double>(input["HeatSourceSpeed"]));
  mhs->initialPosition = CreateEigenVector(py::array_t<double>(input["initialPosition"]));
  mhs->currentPosition    = mhs->initialPosition;

  // TIME DEPENDENCY
  if (input.contains("steadyState")) {
    isSteady = py::cast<bool>(input["steadyState"]);
  }

  // INITIALIZE UNKNOWN
  Tenv = py::cast<double>(input["environmentTemperature"]);
  if (input.contains("depositionTemperature")) {
    Tdeposition = py::cast<double>(input["depositionTemperature"]);
  } else {
    Tdeposition = Tenv;
  }
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
    for (int ifacet : domain.boundaryFacets.getIndices() ) {
      weakBcFacets[ifacet] = 2;
    }
  }

  // ADVECTION
  if (input.contains("advectionSpeed")) {
    advectionSpeed = CreateEigenVector(py::array_t<double>(input["advectionSpeed"]));
    if (advectionSpeed.norm() > 1e-10) isAdvection = true;
  }

  // DOMAIN MOTION
  // TODO: Move this down to mesh level!
  if (input.contains("speedFRF")) {
    domain.mesh->speedFRF = CreateEigenVector(py::array_t<double>(input["speedFRF"]));
  }

  // ASSS STABILIZATION
  if (input.contains("isStabilized")) {
    isStabilized = py::cast<bool>(input["isStabilized"]);
  }
}

Eigen::VectorXd CreateEigenVector(py::array_t<double> n) {
  // Convert 1d numpy vector to 1d eigen vector
  py::buffer_info buf = n.request();
  double* ptr = (double*)buf.ptr;
  Eigen::VectorXd e(buf.size);
  for (size_t idx = 0; idx < buf.size; ++idx) {
    e(idx) = ptr[idx];
  }
  return e;
}
