#include <iostream>
#include <vector>
#include "includes/line.h"
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
  double m_ij, k_ij, ip;

  if (!isAssembled) {
    // initialize data structures
    M.resize(mesh.nnodes, mesh.nnodes); // mass mat
    K.resize(mesh.nnodes, mesh.nnodes); // stiffness mat
    lhs.resize( mesh.nnodes, mesh.nnodes );
    rhs.resize( mesh.nnodes );
    pulse.resize( mesh.nnodes ); // source term

    vector<T> M_coeffs;
    M_coeffs.reserve( 3*mesh.nnodes );
    vector<T> K_coeffs;
    K_coeffs.reserve( 3*mesh.nnodes );

    Line l;
    // matrices assembly
    M.setZero();
    K.setZero();

    for (int ielem = 0; ielem < mesh.nels; ielem++ ) {
      l = mesh.getElement( ielem );
      for (int inode = 0; inode < l.nnodes; inode++) {
        for (int jnode = 0; jnode < l.nnodes; jnode++) {
          m_ij = 0;
          k_ij = 0;
          for (int igp = 0; igp < l.nnodes; igp++) {
            // mass matrix
            m_ij += l.gpweight[igp] * l.baseFunGpVals[inode][igp]*l.baseFunGpVals[jnode][igp]*l.vol;
            // stiffness matrix
            ip = inner_product(l.baseFunGradGpVals[inode][igp].begin(),
                  l.baseFunGradGpVals[inode][igp].end(),
                  l.baseFunGradGpVals[jnode][igp].begin(),
                  0.0);
            k_ij += l.gpweight[igp] * ip * l.vol;
          }
          m_ij *= material["rho"]*material["cp"];
          k_ij *= material["k"];

          // lookup i or j belong to fixed nodes
          //
          // if belong
          // send it to rhs multiplying dirichlet


          M_coeffs.push_back( T( l.con[inode], l.con[jnode], m_ij ) );
          K_coeffs.push_back( T( l.con[inode], l.con[jnode], k_ij ) );
        }
      }
    }
    M.setFromTriplets( M_coeffs.begin(), M_coeffs.end() );
    K.setFromTriplets( K_coeffs.begin(), K_coeffs.end() );
    isAssembled = true;
  }

  //SOLVE
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

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

  lhs.setZero();
  rhs.setZero();
  mhs.computePulse(pulse, time+dt, mesh);

  // general treatment implicit schemes
  lhs += K;
  rhs += pulse;
  switch (timeIntegrator.currentIntegrator) {
    case 0:
      {//Forward Euler
       //Special treatment explicit scheme
       //Undo computations for implicit scheme business
        lhs.setZero();
        rhs.setZero();
        pulse.setZero();
        lhs += M;
        mhs.computePulse(pulse, time, mesh);
        rhs += pulse - K*solution;
        solver.compute(lhs);//overkill
        solution += dt * solver.solve(rhs);
        break;
      }
    default:
      {//Implicit scheme: Generalized treatment

        lhs += timeIntegrator.lhsCoeff * M / dt;
        rhs += M * (prevSolutions(Eigen::placeholders::all, Eigen::seq( 0, timeIntegrator.rhsCoeff.size() - 1)) * timeIntegrator.rhsCoeff) / dt;
        solver.compute( lhs );
        solution = solver.solve(rhs);
        break;
      }
  }

  // End iteration operations
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
        .def_readonly("mhs", &Problem::mhs)
        .def_readonly("mesh", &Problem::mesh)
        .def_readwrite("time", &Problem::time)
        .def_readonly("dt", &Problem::dt);
    py::class_<Mesh>(m, "Mesh")
        .def(py::init<>())
        .def_readonly("pos", &Mesh::pos)
        .def_readonly("nels", &Mesh::nels)
        .def_readonly("nnodes", &Mesh::nnodes)
        .def("initialize1DMesh", &Mesh::initialize1DMesh);
    py::class_<HeatSource>(m, "HeatSource")
        .def(py::init<>())
        .def_readwrite("currentPosition", &HeatSource::currentPosition);
}
