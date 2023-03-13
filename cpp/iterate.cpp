#include <iostream>
#include "Problem.h"

void Problem::iterate() {
  // BEGIN ITERATION
  preIterate();
  // ASSEMBLY
  // LHS, space
  assembleSpatialLHS();
  // RHS, space
  assembleSpatialRHS();
  // LHS & RHS, VMS
  if (isStabilized) {
    assembleStabilization();
  }

  // LHS & RHS, time
  assembleTime();

  // LHS & RHS, inactive nodes
  forceInactiveNodes();

  //Dirichlet BC
  forceDirichletNodes();

  //SOLVE
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  lhs.setFromTriplets( lhsCoeffs.begin(), lhsCoeffs.end() );
  //Solve linear system
  solver.compute( lhs );
  if (not(solver.info() == Eigen::Success)) {
    std::cout << "Singular matrix!" << std::endl;
  }
  unknown.values = solver.solve(rhs);

  //END ITERATION
  postIterate();
}
