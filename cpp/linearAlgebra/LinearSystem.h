#ifndef LINEARSYSTEM
#define LINEARSYSTEM
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <iostream>

// Attempt at distributed linear system
class LinearSystem {
  public:
    Eigen::SparseMatrix<double> lhs;
    Eigen::VectorXd rhs;
    int _ndofs;
    int _rank = 0;
    std::vector<int> dofNumbering;
    std::vector<Eigen::Triplet<double>> lhsCoeffs;
    Eigen::VectorXd sol;

    LinearSystem() = default;
    LinearSystem(int ndofs ) {
      _ndofs = ndofs;
      dofNumbering.resize(_ndofs);
      for (int i = 0; i < _ndofs; ++i) {
        dofNumbering[i] = i;
      }
      _rank = 0;
      if (_rank == 0) {
        allocate();
      }
    }

    void assemble() {
      lhs.setFromTriplets( lhsCoeffs.begin(), lhsCoeffs.end() );
    }
    void cleanup() {
      lhs.setZero();
      rhs.setZero();
      lhsCoeffs.clear();
    }
    void solve() {
      Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
      //Solve linear system
      solver.compute( lhs );
      if (not(solver.info() == Eigen::Success)) {
        std::cout << "Singular matrix!" << std::endl;
        exit(-1);
      }
      sol = solver.solve(rhs);
    }
    int getNdofs() { return _ndofs; }
    int getRank() { return _rank; }
  private:
    void allocate() {
      lhs.resize( _ndofs, _ndofs );
      rhs.resize( _ndofs );
    }
};
#endif
