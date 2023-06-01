#ifndef LINEARSYSTEM
#define LINEARSYSTEM
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <iostream>

class Problem;// Forward declaration of problem class
              //
class LinearSystem {
  public:
    Eigen::SparseMatrix<double> lhs;
    Eigen::VectorXd rhs;
    int _ndofs;
    std::vector<Eigen::Triplet<double>> lhsCoeffs;
    Eigen::VectorXd sol;

    LinearSystem() = default;
    LinearSystem(Problem &p);
    LinearSystem(Problem &p1, Problem &p2);

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
};
#endif
