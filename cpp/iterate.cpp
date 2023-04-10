#include <iostream>
#include "Problem.h"

void Problem::iterate() {
  // BEGIN ITERATION
  preIterate();

  // ASSEMBLY
  // LHS & RHS, inactive nodes
  forceInactiveNodes();
  // LHS, space
  assembleSpatialLHS();
  // RHS, space
  assembleSpatialRHS();

  // Neumann BC
  assembleNeumann();

  // Convection BC
  assembleConvectionLHS();
  assembleConvectionRHS();

  // LHS & RHS, VMS
  if (isStabilized) {
    assembleStabilization();
  }

  // LHS & RHS, time
  assembleTime();

  // Proper assembly
  lhs.setFromTriplets( lhsCoeffs.begin(), lhsCoeffs.end() );

  // Dirichlet BC
  forceDirichletNodes();

  //SOLVE
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  //Solve linear system
  solver.compute( lhs );
  if (not(solver.info() == Eigen::Success)) {
    std::cout << "Singular matrix!" << std::endl;
    exit(-1);
  }
  unknown.values = solver.solve(rhs);

  //END ITERATION
  postIterate();
}
