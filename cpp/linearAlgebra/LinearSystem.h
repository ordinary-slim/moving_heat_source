#ifndef LINEARSYSTEM
#define LINEARSYSTEM
#include "Solver.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <iostream>
#include <memory>

class Problem;// Forward declaration of problem class
namespace mesh { class Domain; }// Forward declaration within a namespace

class LinearSystem 
  : public std::enable_shared_from_this<LinearSystem>
{
  public:
    Eigen::SparseMatrix<double, Eigen::RowMajor> lhs;
    Eigen::VectorXd rhs;
    size_t _ndofs;
    std::vector<Eigen::Triplet<double>> lhsCoeffs;
    Eigen::VectorXd sol;
    Eigen::VectorXd initialGuess;

    LinearSystem() = default;
    LinearSystem(Problem &p);
    LinearSystem(Problem &p1, Problem &p2);

    // For solving projections
    LinearSystem(mesh::Domain &d);

    static std::shared_ptr<LinearSystem> Create(Problem &p1, Problem &p2);

    void cleanup() {
      lhs.setZero();
      rhs.setZero();
      lhsCoeffs.clear();
    }

    void concatenateProblem(Problem &p);
    void concatenateDomain(mesh::Domain &d);

    int getNdofs() { return _ndofs; }
    void allocate() {
      lhs.resize( _ndofs, _ndofs );
      rhs.resize( _ndofs );
    }

    void setInitialGuess(Problem* p1, Problem* p2 = nullptr);
    void assemble();

    std::unique_ptr<Solver> solver = std::make_unique<EigenBiCGSTAB_IncompleteLUT>();
    void solve();
    void setSolver(int idxSolver);

};
#endif
