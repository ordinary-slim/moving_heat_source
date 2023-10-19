#ifndef LINEARSYSTEM
#define LINEARSYSTEM
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <iostream>
#include <memory>

class Problem;// Forward declaration of problem class
              //
void solveEigenBiCGSTAB( Eigen::SparseMatrix<double> &lhs,
                         Eigen::VectorXd &rhs,
                         Eigen::VectorXd &sol );
void solveEigenCG( Eigen::SparseMatrix<double> &lhs,
                   Eigen::VectorXd &rhs,
                   Eigen::VectorXd &sol );


class LinearSystem 
  : public std::enable_shared_from_this<LinearSystem>
{
  public:
    Eigen::SparseMatrix<double> lhs;
    Eigen::VectorXd rhs;
    size_t _ndofs;
    std::vector<Eigen::Triplet<double>> lhsCoeffs;
    Eigen::VectorXd sol;

    LinearSystem() = default;
    LinearSystem(Problem &p);
    LinearSystem(Problem &p1, Problem &p2);

    static std::shared_ptr<LinearSystem> Create(Problem &p1, Problem &p2);

    void cleanup() {
      lhs.setZero();
      rhs.setZero();
      lhsCoeffs.clear();
    }

    void concatenateProblem(Problem &p);

    int getNdofs() { return _ndofs; }
    void allocate() {
      lhs.resize( _ndofs, _ndofs );
      rhs.resize( _ndofs );
    }

    void assemble();
    void solve();
    void (*externalSolve)( Eigen::SparseMatrix<double> &,
                           Eigen::VectorXd &,
                           Eigen::VectorXd &) = &solveEigenBiCGSTAB;
    void setSolver(bool isSymmetric = false);
};
#endif
