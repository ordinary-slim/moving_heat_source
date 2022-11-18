#include <iostream>
#include <vector>
#include "line.h"
#include "problem.h"
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

  // initialize data structures
  SpMat M(mesh.nnodes, mesh.nnodes); // mass mat
  vector<T> M_coeffs;
  M_coeffs.reserve( 3*mesh.nnodes );
  SpMat K(mesh.nnodes, mesh.nnodes); // stiffness mat
  vector<T> K_coeffs;
  K_coeffs.reserve( 3*mesh.nnodes );
  Eigen::VectorXd pulse( mesh.nnodes ); // source term

  Line l;
  // matrices assembly
  M.setZero();
  K.setZero();

  vector<int> boundaryNodes;

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


  //slice
  //node 0 and node nnodes-1 are dirichlet nodes at IC

  //solve
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

  //load vector computation/assembly inside of time integration
  pulse.setZero();

  if (timeIntegration == "ForwardEuler") {
    solver.compute(M);//overkill
    Eigen::VectorXd diffusion = K*solution;
    mhs.computePulse(pulse, time, mesh);
    Eigen::VectorXd rhs = pulse - diffusion;
    deltaSolution = dt * solver.solve(rhs);

    solution += deltaSolution;
  } else if (timeIntegration == "BackwardEuler") {
    solver.compute(M + dt * K);
    Eigen::VectorXd diffusion = K*solution;
    Eigen::VectorXd rhs = M * solution + dt * pulse; 
    deltaSolution = dt * solver.solve(rhs);

    solution += deltaSolution;
  } else {
    std::cout << "Check time integration string" << std::endl;
    std::exit(1);
  }

  time += dt;
/*
if (plot_source) {
  //convert from Eigen::Vector to std::vector
  vector<double> vpulse( pulse.data(), pulse.data() + pulse.size());

  // get power density
  vector<double> pd( fineMesh.pos.size() );
  transform( fineMesh.pos.begin(), fineMesh.pos.end(), pd.begin(),
      [&x0, &radius, &P](double x) { return power_density( x0, x, radius, P ); } );

  plt::clf();
  plt::plot( std::vector<double>({x0}), std::vector<double>({0.0}),
      {{"linestyle", "none"}, {"marker", "X"}, {"color", "red"}, {"markersize", "20"}, {"label", "x0"}});
  plt::plot( mesh.pos, vpulse, {{"marker", "o"}, {"color", "red"}, {"label", "Pi"}} );
  plt::plot( fineMesh.pos, pd, {{"linestyle", "--"}, {"color", "blue"}, {"label", "f"}} );
  plt::xlim( 0.0, L );
  double maxPower = max( maxPower, *max_element( pulse.begin(), pulse.end() ) );
  //plt::ylim( 0.0, 1.1*maxPower );

  plt::legend();
  plt::pause( 0.04 );
  plt::save( "tmp.png", 300 );
}
if (plot_solution) {
  vector<double> vsolution( solution.data(), solution.data() + solution.size());

  plt::clf();
  plt::plot( std::vector<double>({x0}), std::vector<double>({0.0}),
      {{"linestyle", "none"}, {"marker", "X"}, {"color", "red"}, {"markersize", "20"}, {"label", "x0"}});
  plt::plot( mesh.pos, vsolution, {{"color", "red"}, {"label", "sol"}});
  plt::xlim( 0.0, L );
  plt::ylim( -1.0, 50.0 );

  plt::annotate("t = " + to_string( time ), 0.05, 0.95);
  plt::legend();
  plt::pause( 0.16 );
  plt::save("iter" + to_string(iter) + ".png", 300 );
}
*/

}

PYBIND11_MODULE(MovingHeatSource, m) {
    py::class_<Problem>(m, "Problem", py::dynamic_attr())
        .def(py::init<>())
        .def("initialize", &Problem::initialize)
        .def("iterate", &Problem::iterate)
        .def_readonly("solution", &Problem::solution)
        .def_readonly("mesh", &Problem::mesh)
        .def_readonly("time", &Problem::time)
        .def_readonly("dt", &Problem::dt);
    py::class_<Mesh>(m, "Mesh")
        .def(py::init<>())
        .def_readonly("pos", &Mesh::pos)
        .def("initialize1DMesh", &Mesh::initialize1DMesh);
}
