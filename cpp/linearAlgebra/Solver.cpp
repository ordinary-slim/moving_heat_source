#include "LinearSystem.h"
#include "Solver.h"
#include <unsupported/Eigen/SparseExtra>

template<typename EigenSolver>
void eigenIterativeSolve( EigenSolver &eigenSolver,
    Eigen::SparseMatrix<double, Eigen::RowMajor> &lhs, Eigen::VectorXd &rhs, Eigen::VectorXd &sol, Eigen::VectorXd &initialGuess) {
  eigenSolver.compute( lhs );
  if (not(eigenSolver.info() == Eigen::Success)) {
    std::cout << "Singular matrix!" << std::endl;
    exit(-1);
  }
  if (initialGuess.size()) {
    sol = eigenSolver.solveWithGuess(rhs, initialGuess);
  } else {
    sol = eigenSolver.solve(rhs);
  }
  std::cout << "#Iterations:     " << eigenSolver.iterations() << std::endl;
  std::cout << "Estimated error: " << eigenSolver.error()      << std::endl;

  if (std::isnan(eigenSolver.error())) {
    Eigen::saveMarket(lhs, "faulty_lhs.mtx");
    Eigen::saveMarketVector(rhs, "faulty_rhs.mtx");
    throw std::runtime_error("Solver did not converge, wrote mtx files.");
  }
}

void EigenCG::solve(LinearSystem *ls) {

  Eigen::ConjugateGradient<Eigen::SparseMatrix<double, Eigen::RowMajor> > eigenSolver;
  std::cout << "Using CG" << std::endl;
  eigenIterativeSolve( eigenSolver, ls->lhs, ls->rhs, ls->sol, ls->initialGuess );
}

void EigenBiCGSTAB_Jacobi::solve(LinearSystem *ls) {

  Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>> eigenSolver;
  std::cout << "Using BiCGSTAB with diagonal preconditioner" << std::endl;
  eigenIterativeSolve( eigenSolver, ls->lhs, ls->rhs, ls->sol, ls->initialGuess );
}

void EigenBiCGSTAB_IncompleteLUT::solve(LinearSystem *ls) {

  Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>, Eigen::IncompleteLUT<double> > eigenSolver;
  std::cout << "Using BiCGSTAB with IncompleteLUT preconditioner" << std::endl;
  eigenIterativeSolve( eigenSolver, ls->lhs, ls->rhs, ls->sol, ls->initialGuess );
}
