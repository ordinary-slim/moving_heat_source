#include "includes/problem.h"
#include <map>
#include <string>
#include "../external/pybind11/include/pybind11/pybind11.h"
#include "../external/pybind11/include/pybind11/eigen.h"
#include "../external/pybind11/include/pybind11/stl.h"

namespace py = pybind11;

void Problem::initialize(py::dict &input) {
  // tstepping
  dt = py::cast<double>(input["dt"]);
  // MESH
  // 1D mesh generated by software
  if (input.contains("1D")) {
    double a = py::cast<double>(input["Left"]);
    double b = py::cast<double>(input["Right"]);
    int nels = py::cast<int>( input["nels"] );
    mesh.generate1DMesh( a, b, nels );
  } else {
    //READ POINTS
    py::array points = input["points"];
    int aux_npoints = points.shape(0);
    int aux_ndims =    points.shape(1);
    mesh.pos.resize( aux_npoints, 3 );
    mesh.pos.setZero();
    auto aux_points = points.unchecked<double>();
    for ( int ipoint = 0; ipoint < aux_npoints; ipoint++) {
      for ( int idim = 0; idim < aux_ndims; idim++) {
        mesh.pos(ipoint, idim) =  aux_points(ipoint, idim);
      }
    }
    //READ ELEMENTS
    py::array cells = input["cells"];
    int aux_ncells =       cells.shape(0);
    int aux_nodesPerCell = cells.shape(1);
    mesh.con.resize( aux_ncells, aux_nodesPerCell );
    mesh.con.setOnes();
    mesh.con *= -1;
    auto aux_cells = cells.unchecked<int>();
    for ( int icell = 0; icell < aux_ncells; icell++) {
      for ( int inode = 0; inode < aux_nodesPerCell; inode++) {
        mesh.con(icell, inode) =  aux_cells(icell, inode);
      }
    }
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
