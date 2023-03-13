#include "../Problem.h"
#include <Eigen/Core>
#include <Eigen/Sparse>

void Problem::forceDirichletNodes() {
  // Has to be last step. Similar treatment to Dirichlet BC
  int inode;
  double val;
  for (int iDirichletNode = 0; iDirichletNode < unknown.dirichletNodes.size(); ++iDirichletNode) {
    inode = unknown.dirichletNodes[iDirichletNode];
    val = unknown.dirichletValues[iDirichletNode];
    lhs.coeffRef( inode, inode ) = 1.0;
    rhs[inode] = val;
  }
  lhs.prune( [this](int i, int j, double) {
      return ( (find(unknown.dirichletNodes.begin(), unknown.dirichletNodes.end(), i) == unknown.dirichletNodes.end()) || (i==j) );
      });
}
