#include "../Problem.h"
#include <Eigen/Core>
#include <Eigen/Sparse>

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;

void Problem::forceInactiveNodes() {
  // Similar treatment to Dirichlet BC
  if (domain.hasInactive) {
    // Treat inactive nodes
    vector<T> InacNodes_coeffs;
    InacNodes_coeffs.reserve( domain.mesh->nnodes );
    for (int inode = 0; inode < domain.mesh->nnodes; inode++) {
      if (!domain.activeNodes[inode]) {
        InacNodes_coeffs.push_back( T(inode, inode, 1) );
        rhs[ inode ] = unknown.values[ inode ];
      }
    }
    lhsCoeffs.insert( lhsCoeffs.end(), InacNodes_coeffs.begin(), InacNodes_coeffs.end() );
    massCoeffs.insert( massCoeffs.end(), InacNodes_coeffs.begin(), InacNodes_coeffs.end() );
  }
}
