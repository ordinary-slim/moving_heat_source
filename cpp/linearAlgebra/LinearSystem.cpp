#include "../Problem.h"
#include "../mesh/Domain.h"
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

void LinearSystem::concatenateDomain(mesh::Domain &d) {
  /*
   * Add domain to numbering
   */
  d.dofNumbering.clear();
  d.dofNumbering.reserve(d.mesh->nnodes);

  for (int inode = 0; inode < d.mesh->nnodes; inode++){
    if (not(d.activeNodes[inode])) {
      d.dofNumbering.push_back( -1 );
    } else {
      d.dofNumbering.push_back( _ndofs );
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

LinearSystem::LinearSystem(mesh::Domain &d) {
  _ndofs = 0;
  concatenateDomain( d );
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
  p1.assembling2external = true;
  p2.assembling2external = true;
  return p1.ls;
}

void LinearSystem::solve() {
  if (not(_ndofs)) { return; }
  externalSolve( lhs, rhs, sol, initialGuess );
}

void LinearSystem::assemble() {
  lhs.setFromTriplets( lhsCoeffs.begin(), lhsCoeffs.end() );
}

void LinearSystem::setInitialGuess(Problem* p1, Problem* p2) {
  initialGuess.resize( _ndofs );
  p1->setInitialGuess();
  if (p2 != nullptr) {
    p2->setInitialGuess();
  }
}

void LinearSystem::setSolver(bool isSymmetric) {
  if (isSymmetric) {
    externalSolve = &solveEigenCG;
  } else {
    externalSolve = &solveEigenBiCGSTAB;
  }
}

void solveEigenBiCGSTAB( Eigen::SparseMatrix<double, Eigen::RowMajor> &lhs,
                         Eigen::VectorXd &rhs,
                         Eigen::VectorXd &sol,
                         Eigen::VectorXd &initialGuess) {

  Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor> > solver;
  solver.compute( lhs );
  if (not(solver.info() == Eigen::Success)) {
    std::cout << "Singular matrix!" << std::endl;
    exit(-1);
  }
  if (initialGuess.size()) {
    std::cout << "Solving with guess " << std::endl;
    sol = solver.solveWithGuess(rhs, initialGuess);
  } else {
    sol = solver.solve(rhs);
  }
  std::cout << "EigenBiCGSTAB #iterations:     " << solver.iterations() << std::endl;
  std::cout << "EigenBiCGSTAB estimated error: " << solver.error()      << std::endl;
}
void solveEigenCG( Eigen::SparseMatrix<double, Eigen::RowMajor> &lhs,
                   Eigen::VectorXd &rhs,
                   Eigen::VectorXd &sol,
                   Eigen::VectorXd &initialGuess) {

  Eigen::ConjugateGradient<Eigen::SparseMatrix<double, Eigen::RowMajor> > solver;
  solver.compute( lhs );
  if (not(solver.info() == Eigen::Success)) {
    std::cout << "Singular matrix!" << std::endl;
    exit(-1);
  }
  if (initialGuess.size()) {
    sol = solver.solveWithGuess(rhs, initialGuess);
  } else {
    sol = solver.solve(rhs);
  }
  std::cout << "CG #iterations:     " << solver.iterations() << std::endl;
  std::cout << "CG estimated error: " << solver.error()      << std::endl;
}
