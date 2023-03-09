#include <iostream>
#include <vector>
#include "mesh/Element.h"
#include "Problem.h"
#include <numeric>
#include <algorithm>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <string>

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


  //SOLVE
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  //Solve linear system
  solver.compute( lhs );
  if (not(solver.info() == Eigen::Success)) {
    std::cout << "Singular matrix!" << std::endl;
  }
  unknown.values = solver.solve(rhs);

  //END ITERATION
  postIterate();
}
