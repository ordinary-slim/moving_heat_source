#include <iostream>
#include "Problem.h"

void Problem::solve() {
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  //Solve linear system
  solver.compute( lhs );
  if (not(solver.info() == Eigen::Success)) {
    std::cout << "Singular matrix!" << std::endl;
    exit(-1);
  }
  unknown.values = solver.solve(rhs);
}
