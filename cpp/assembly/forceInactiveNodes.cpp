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
        ls->lhsCoeffs.push_back( T(dofNumbering[inode], dofNumbering[inode], 1) );
        domain.massCoeffs.push_back( T(inode, inode, 1) );
        ls->rhs[ dofNumbering[inode] ] = unknown.values[ inode ];
      }
    }
  }
}
