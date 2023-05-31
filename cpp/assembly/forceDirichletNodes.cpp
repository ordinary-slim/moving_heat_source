#include "../Problem.h"
#include <Eigen/Core>
#include <Eigen/Sparse>

void Problem::forceDirichletNodes() {
  // Has to be last step. Similar treatment to Dirichlet BC
  // Does it have to be last step?
  int inode;
  double val;
  vector<int> dirichletIndices = dirichletNodes.getTrueIndices();
  for (int inode : dirichletIndices) {
    val = dirichletValues[inode];
    ls.lhs.coeffRef( ls.dofNumbering[inode], ls.dofNumbering[inode] ) = 1.0;
    ls.rhs[ls.dofNumbering[inode]] = val;
  }
  vector<int> globalDirichletIndices( dirichletIndices.size() );
  for (int idx = 0; idx < dirichletIndices.size(); ++idx) {
    globalDirichletIndices[idx] = ls.dofNumbering[ dirichletIndices[idx] ];
  }
  ls.lhs.prune( [&globalDirichletIndices](int i, int j, double) {
      return ( (find(globalDirichletIndices.begin(), globalDirichletIndices.end(), i) == globalDirichletIndices.end()) || (i==j) );
      });
}
