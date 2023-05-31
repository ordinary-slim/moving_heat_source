#include "../Problem.h"
#include "LinearSystem.h"

LinearSystem::LinearSystem(Problem &p) {
  _ndofs = p.domain.mesh->nnodes;
  p.dofNumbering.resize(_ndofs);
  for (int i = 0; i < _ndofs; ++i) {
    p.dofNumbering[i] = i;
  }
  allocate();
}
LinearSystem::LinearSystem(Problem &p1, Problem &p2) {
  /*
   * Linear system shared by two problems.
   * Could be a vector of problems but don't know how
   * to pass by reference from Python side
   */
  _ndofs = p1.domain.mesh->nnodes + p2.domain.mesh->nnodes ;
  p1.dofNumbering.resize( p1.domain.mesh->nnodes );
  for (int inode = 0; inode < p1.domain.mesh->nnodes; ++inode ) {
    p1.dofNumbering[inode] = inode;
  }
  p1.ls = this;
  p1.assembling2external = true;
  p2.dofNumbering.resize( p2.domain.mesh->nnodes );
  for (int inode = p1.domain.mesh->nnodes; inode < _ndofs; ++inode ) {
    p2.dofNumbering[inode-p1.domain.mesh->nnodes] = inode;
  }
  p2.ls = this;
  p2.assembling2external = true;
  allocate();
}
void LinearSystem::solve() {
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  //Solve linear system
  solver.compute( lhs );
  if (not(solver.info() == Eigen::Success)) {
    std::cout << "Singular matrix!" << std::endl;
    exit(-1);
  }
  sol = solver.solve(rhs);
}

void LinearSystem::assemble() {
  lhs.setFromTriplets( lhsCoeffs.begin(), lhsCoeffs.end() );
  forceDofs();
}

void LinearSystem::forceDofs() {
  // Has to be last step. Similar treatment to Dirichlet BC
  // Does it have to be last step?
  for (int i = 0; i < indicesForcedDofs.size(); ++i) {
    int inode = indicesForcedDofs[i];
    double val = valuesForcedDofs[i];
    lhs.coeffRef( inode, inode ) = 1.0;
    rhs[inode] = val;
  }
  lhs.prune( [this](int i, int j, double) {
      return ( (find(indicesForcedDofs.begin(), indicesForcedDofs.end(), i) == indicesForcedDofs.end()) || (i==j) );
      });
}
