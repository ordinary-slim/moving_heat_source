#include "../Problem.h"
#include "LinearSystem.h"

LinearSystem::LinearSystem(Problem &p) {
  p.dofNumbering.clear();
  p.dofNumbering.reserve(p.domain.mesh->nnodes);
  _ndofs = 0;
  for (int inode = 0; inode < p.domain.mesh->nnodes; inode++){
    if (p.forcedDofs[inode] == 0 ) {
      p.dofNumbering.push_back( _ndofs );
      ++_ndofs;
    } else {
      p.dofNumbering.push_back( -1 );
    }
  }
  allocate();
  cleanup();
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
}
