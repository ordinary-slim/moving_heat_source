#include <iostream>
#include <vector>
#include "includes/element.h"
#include "includes/problem.h"
#include <numeric>
#include <algorithm>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <string>
#include "../external/pybind11/include/pybind11/pybind11.h"
#include "../external/pybind11/include/pybind11/stl.h"
#include "../external/pybind11/include/pybind11/eigen.h"

using namespace std;
namespace py = pybind11;

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;

void Problem::iterate() {
  // numerical params
  double m_ij, k_ij, a_ij, ip;

  if (!isAssembled) {
    // initialize data structures
    M.resize(mesh.nnodes, mesh.nnodes); // mass mat
    K.resize(mesh.nnodes, mesh.nnodes); // stiffness mat
    A.resize(mesh.nnodes, mesh.nnodes); // advection mat
    lhs.resize( mesh.nnodes, mesh.nnodes );
    rhs.resize( mesh.nnodes );
    pulse.resize( mesh.nnodes ); // source term

    vector<T> M_coeffs;
    M_coeffs.reserve( 3*mesh.nnodes );
    vector<T> K_coeffs;
    K_coeffs.reserve( 3*mesh.nnodes );
    vector<T> A_coeffs;
    A_coeffs.reserve( 3*mesh.nnodes );

    // matrices assembly
    M.setZero();
    K.setZero();
    A.setZero();

    Element e;
    double rho = material["rho"];
    double cp = material["cp"];
    double k = material["k"];
    for (int ielem = 0; ielem < mesh.nels; ielem++ ) {
      e = mesh.getElement( ielem );
      for (int inode = 0; inode < e.nnodes; inode++) {
        for (int jnode = 0; jnode < e.nnodes; jnode++) {
          m_ij = 0;
          k_ij = 0;
          a_ij = 0;
          for (int igp = 0; igp < e.nnodes; igp++) {
            // mass matrix
            m_ij += e.gpweight[igp] * e.BaseGpVals[inode][igp]*e.BaseGpVals[jnode][igp]*e.vol;
            // stiffness matrix
            ip = e.GradBaseGpVals[inode][igp].dot( e.GradBaseGpVals[jnode][igp] );
            k_ij += e.gpweight[igp] * ip * e.vol;
            // advection matrix
            ip = e.GradBaseGpVals[jnode][igp].dot(advectionSpeed);
            a_ij += e.gpweight[igp] * (ip * e.BaseGpVals[inode][igp]) * e.vol;
          }
          m_ij *= rho * cp ;
          k_ij *= k ;
          a_ij *= rho * cp;

          // lookup i or j belong to fixed nodes
          //
          // if belong
          // send it to rhs multiplying dirichlet


          M_coeffs.push_back( T( e.con[inode], e.con[jnode], m_ij ) );
          K_coeffs.push_back( T( e.con[inode], e.con[jnode], k_ij ) );
          A_coeffs.push_back( T( e.con[inode], e.con[jnode], a_ij ) );
        }
      }
    }
    M.setFromTriplets( M_coeffs.begin(), M_coeffs.end() );
    K.setFromTriplets( K_coeffs.begin(), K_coeffs.end() );
    A.setFromTriplets( A_coeffs.begin(), A_coeffs.end() );
    isAssembled = true;
  }

  //SOLVE
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  lhs.setZero();
  rhs.setZero();

  // general treatment implicit schemes
  mhs.computePulse(pulse, mesh, time+dt, dt);
  lhs += K;
  if (isAdvection) lhs += A;
  rhs += pulse;

  if (not isSteady) {
    //Set time integration
    if (timeIntegrator.nstepsStored < timeIntegrator.nstepsRequired ) {
      if (timeIntegrator.nstepsStored >= 4) {
        timeIntegrator.setCurrentIntegrator( 4 );
      } else if (timeIntegrator.nstepsStored >= 3) {
        timeIntegrator.setCurrentIntegrator( 3 );
      } else if (timeIntegrator.nstepsStored >= 2) {
        timeIntegrator.setCurrentIntegrator( 2 );
      } else {
        timeIntegrator.setCurrentIntegrator( 1 );
      }
    } else {
      timeIntegrator.setCurrentIntegrator( timeIntegrator.desiredIntegrator );
    }
    //Add time dependency
    lhs += timeIntegrator.lhsCoeff * M / dt;
    rhs += M * (prevSolutions(Eigen::placeholders::all, Eigen::seq( 0, timeIntegrator.rhsCoeff.size() - 1)) * timeIntegrator.rhsCoeff) / dt;
  }

  //Solve linear system
  solver.compute( lhs );
  if (not(solver.info() == Eigen::Success)) {
    cout << "Singular matrix!" << endl;
  }
  solution = solver.solve(rhs);

  //END ITERATION
  postIterate();

}

PYBIND11_MODULE(MovingHeatSource, m) {
    py::class_<Problem>(m, "Problem", py::dynamic_attr())
        .def(py::init<>())
        .def("initialize", &Problem::initialize)
        .def("initializeIntegrator", &Problem::initializeIntegrator)
        .def("iterate", &Problem::iterate)
        .def("setTime", &Problem::setTime)
        .def_readwrite("solution", &Problem::solution)
        .def_readwrite("pulse", &Problem::pulse)//debugging
        .def_readonly("mhs", &Problem::mhs)
        .def_readonly("mesh", &Problem::mesh)
        .def_readwrite("time", &Problem::time)
        .def_readonly("dt", &Problem::dt)
        .def_readonly("isAdvection", &Problem::isAdvection)
        .def_readonly("advectionSpeed", &Problem::advectionSpeed);
    py::class_<Mesh>(m, "Mesh", py::dynamic_attr())
        .def(py::init<>())
        .def_readonly("pos", &Mesh::pos)
        .def_readonly("nels", &Mesh::nels)
        .def_readonly("nnodes", &Mesh::nnodes)
        .def("generate1DMesh", &Mesh::generate1DMesh)
        .def("getElement", &Mesh::getElement);
    py::class_<Element>(m, "Element", py::dynamic_attr())
        .def(py::init<>())
        .def_readonly("pos", &Element::pos)
        .def_readonly("nnodes", &Element::nnodes)
        .def_readonly("gpos", &Element::gpos)
        .def_readonly("ngpoints", &Element::ngpoints)
        .def_readonly("gpweight", &Element::gpweight)
        .def_readonly("con", &Element::con)
        .def_readonly("vol", &Element::vol)
        .def_readonly("dimension", &Element::dim)
        .def_readonly("elementType", &Element::elementType)
        .def_readonly("GradBaseGpVals", &Element::GradBaseGpVals)
        .def("computeDerivatives", &Element::computeDerivatives);
    py::class_<HeatSource>(m, "HeatSource")
        .def(py::init<>())
        .def_readwrite("currentPosition", &HeatSource::currentPosition);
}
