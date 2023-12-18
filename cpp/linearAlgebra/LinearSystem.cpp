#include "../Problem.h"
#include "../mesh/Domain.h"
#include "LinearSystem.h"
#include "Solver.h"
#include "MATLAB.h"
#include <cmath>

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
  solver->solve( this );
}

void LinearSystem::assemble() {
  lhs.setFromTriplets( lhsCoeffs.begin(), lhsCoeffs.end() );
  lhsCoeffs.clear();
}

void LinearSystem::setInitialGuess(Problem* p1, Problem* p2) {
  initialGuess.resize( _ndofs );
  p1->setInitialGuess();
  if (p2 != nullptr) {
    p2->setInitialGuess();
  }
}

void LinearSystem::setSolver(int idxSolver) {
  switch (idxSolver) {
    case 1:
      solver = std::make_unique<EigenCG>();
      break;
    case 2:
      solver = std::make_unique<EigenBiCGSTAB_Jacobi>();
      break;
    case 3:
      solver = std::make_unique<MatLabSolver>();
      break;
    case 0: default:
      solver = std::make_unique<EigenBiCGSTAB_IncompleteLUT>();
      break;
  }
}
