#include "../Problem.h"
#include <Eigen/Core>
#include <Eigen/Sparse>

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;

void Problem::forceInactiveNodes() {
  // Has to be last step. Similar treatment to Dirichlet BC
  if (mesh.hasInactive) {
    // Treat inactive nodes
    vector<T> InacNodes_coeffs;
    InacNodes_coeffs.reserve( mesh.nnodes );
    for (int inode = 0; inode < mesh.nnodes; inode++) {
      if (mesh.activeNodes[inode] == 0) {
        InacNodes_coeffs.push_back( T(inode, inode, 1) );
        rhs[ inode ] = unknown.values[ inode ];
      }
    }
    I.setFromTriplets( InacNodes_coeffs.begin(), InacNodes_coeffs.end() );
    lhs += I;
  }
}
