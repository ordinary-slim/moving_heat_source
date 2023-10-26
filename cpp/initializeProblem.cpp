#include "Problem.h"
#include <map>
#include <string>
#include <algorithm>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

namespace py = pybind11;

Problem::Problem(mesh::Mesh &mesh, py::dict &input) :
  domain( mesh::Domain( &mesh , this ) ),
  material( input ),
  unknown( fem::Function( &domain ) ),
  forcedDofs( mesh::MeshTag<int>( &mesh, 0, 0)),
  dirichletNodes( mesh::MeshTag<int>( &mesh ) ),
  dirichletValues( mesh::MeshTag<double>( &mesh ) ),
  weakBcFacets( mesh::MeshTag<int>( &mesh, mesh.dim-1 ) ),
  neumannFluxes( mesh::MeshTag<vector<double>>( &mesh, mesh.dim-1, vector<double>() ) ),
  gammaNodes( mesh::MeshTag<int>( &mesh, 0 ) ),
  gammaFacets( mesh::MeshTag<int>( &mesh, mesh.dim-1 ) )
{

  // HEAT SOURCE
  // set type of source term
  switch (int(py::cast<int>( input["sourceTerm"] ))) {
    case 1:
      mhs = std::make_unique<heat::gaussianPowerDensity1D>(input, this);
      break;
    case 2 :
      mhs = std::make_unique<heat::gaussianPowerDensity2D>(input, this);
      break;
    case 3 :
      mhs = std::make_unique<heat::gaussianPowerDensity3D>(input, this);
      break;
    case 11: { 
        double heatSouceWidth  = py::cast<double>(input["heatSourceWidth"]);
        double heatSouceHeight = py::cast<double>(input["heatSourceHeight"]);
        mhs = std::make_unique<heat::LumpedHeatSource>( heatSouceWidth, heatSouceHeight, input, this);
        break; }
    case 86: {
        mhs = std::make_unique<heat::cteHeat>(input, this);
        break; }
    default:
      {
        if (domain.mesh->dim == 1 ) {
          mhs = std::make_unique<heat::gaussianPowerDensity1D>(input, this);
        } else if (domain.mesh->dim == 2 ) {
          mhs = std::make_unique<heat::gaussianPowerDensity2D>(input, this);
        } else {
          mhs = std::make_unique<heat::gaussianPowerDensity3D>(input, this);
        }
        break;
      }
  }

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
  if (material.convectionCoeff > 0.0) {
    for (int ifacet : domain.boundaryFacets.getIndices() ) {
      weakBcFacets[ifacet] = 2;
    }
  }

  // ADVECTION
  if (input.contains("advectionSpeed")) {
    advectionSpeed = CreateEigenVector(py::array_t<double>(input["advectionSpeed"]));
    if (advectionSpeed.norm() > 1e-10) isAdvection = true;
  }

  // STABILIZATION
  setStabilization( input );

  // DOMAIN MOTION
  // TODO: Move this down to mesh level!
  if (input.contains("speedDomain")) {
    domain.setSpeed( CreateEigenVector(py::array_t<double>(input["speedDomain"])) );
  }

  // Choose symmetric solver if applies
  if (input.contains("isSymmetric")) {
    isSymmetric = py::cast<bool>(input["isSymmetric"]);
  }
}

Eigen::VectorXd CreateEigenVector(py::array_t<double> n) {
  // Where should I put this?
  py::buffer_info buf = n.request();
  double* ptr = (double*)buf.ptr;
  Eigen::VectorXd e(buf.size);
  for (size_t idx = 0; idx < buf.size; ++idx) {
    e(idx) = ptr[idx];
  }
  return e;
}
