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
    lhs.coeffRef( inode, inode ) = 1.0;
    rhs[inode] = val;
  }
  lhs.prune( [&dirichletIndices](int i, int j, double) {
      return ( (find(dirichletIndices.begin(), dirichletIndices.end(), i) == dirichletIndices.end()) || (i==j) );
      });
}
