#include <iostream>
#include <vector>
#include "domain/element.h"
#include "problem.h"
#include <numeric>
#include <algorithm>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <string>

using namespace std;
typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;

void Problem::iterate() {
  // ASSEMBLY
  assemble();

  //SOLVE
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  lhs.setZero();
  rhs.setZero();

  // general treatment implicit schemes
  mhs.computePulse(pulse, mesh, time+dt, dt);//assembly of RHS
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

  if (mesh.hasInactive) {
    // Treat inactive nodes
    vector<T> InacNodes_coeffs;
    InacNodes_coeffs.reserve( mesh.nnodes );
    for (int inode = 0; inode < mesh.nnodes; inode++) {
      if (mesh.activeNodes[inode] == 0) {
        InacNodes_coeffs.push_back( T(inode, inode, 1) );
        rhs[ inode ] = solution[ inode ];
      }
    }
    I.setFromTriplets( InacNodes_coeffs.begin(), InacNodes_coeffs.end() );
    lhs += I;
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
