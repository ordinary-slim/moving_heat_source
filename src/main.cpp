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

double power_density( double x0, double x, double radius, double P );
void compute_pulse(Eigen::VectorXd &pulse, double x0, double radius, double P, Mesh &m);

void Problem::iterate() {
  // params
  double radius=2.0, P=1000000.0;
  double x0 = 5.0, speed=10, Tfinal=10.0;
  // numerical params
  double CFL;
  double m_ij, k_ij, ip;

  // set timestep
  /*
  CFL = 5.0;
  dt = CFL * (L / nels) / speed;
  */
  //maxIter = L / (speed*dt) + 5;

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
        M_coeffs.push_back( T( l.con[inode], l.con[jnode], m_ij ) );
        K_coeffs.push_back( T( l.con[inode], l.con[jnode], k_ij ) );
      }
    }
  }
  M.setFromTriplets( M_coeffs.begin(), M_coeffs.end() );
  K.setFromTriplets( K_coeffs.begin(), K_coeffs.end() );

  //load vector computation/assembly
  pulse.setZero();
  compute_pulse(pulse, x0, radius, P, mesh);

  //solve
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(M);//overkill
  Eigen::VectorXd rhs = pulse + K * solution;
  deltaSolution = dt * solver.solve(rhs);

  solution += deltaSolution;

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

  time += dt;
  //update laser position
  x0 = x0 + dt * speed;
}

PYBIND11_MODULE(MovingHeatSource, m) {
    py::class_<Problem>(m, "Problem", py::dynamic_attr())
        .def(py::init<>())
        .def("initialize", &Problem::initialize)
        .def("iterate", &Problem::iterate)
        .def_readonly("solution", &Problem::solution)
        .def_readonly("mesh", &Problem::mesh)
        .def_readonly("dt", &Problem::time)
        .def_readonly("time", &Problem::dt);
    py::class_<Mesh>(m, "Mesh")
        .def(py::init<>())
        .def_readonly("pos", &Mesh::pos)
        .def("initialize1DMesh", &Mesh::initialize1DMesh);
}
