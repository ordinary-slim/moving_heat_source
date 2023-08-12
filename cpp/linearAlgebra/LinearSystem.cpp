#include "../Problem.h"
#include "LinearSystem.h"

void LinearSystem::concatenateProblem(Problem &p) {
  /*
   * Add problem to numbering
   */
  p.dofNumbering.clear();
  p.dofNumbering.reserve(p.domain.mesh->nnodes);
  p.freeDofsNumbering.clear();
  p.freeDofsNumbering.reserve(p.domain.mesh->nnodes);

  for (int inode = 0; inode < p.domain.mesh->nnodes; inode++){
    if (p.forcedDofs[inode] == 0 ) {
      p.dofNumbering.push_back( _ndofs );
      p.freeDofsNumbering.push_back( _ndofs );
      ++_ndofs;
    } else if (p.forcedDofs[inode] == 1) {
      p.dofNumbering.push_back( -1 );
      p.freeDofsNumbering.push_back( -1 );
    } else if (p.forcedDofs[inode] == 2) {
      // Gamma Dirichlet nodes, special treatment
      p.dofNumbering.push_back( _ndofs );
      p.freeDofsNumbering.push_back( -1 );
      ++_ndofs;
    }
  }
}

LinearSystem::LinearSystem(Problem &p) {
  _ndofs = 0;
  concatenateProblem( p );
  allocate();
  cleanup();
}

LinearSystem::LinearSystem(Problem &p1, Problem &p2) {
  /*
   * Linear system shared by two problems.
   * Could be a vector of problems but don't know how
   * to pass by reference from Python side
   */
  _ndofs = 0 ;

  concatenateProblem( p1 );
  concatenateProblem( p2 );

  p1.assembling2external = true;//this could be removed
  p2.assembling2external = true;//this could be removed

  allocate();
  cleanup();
}

std::shared_ptr<LinearSystem> LinearSystem::Create(Problem &p1, Problem &p2) {
  // Wrapping constructor. Is this good idea?
  p1.ls = std::make_shared<LinearSystem>(p1, p2);
  p2.ls = p1.ls->shared_from_this();

  return p1.ls;
}

void LinearSystem::solve() {
  //Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
  //Solve linear system
  if (not(_ndofs)) { return; }
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
