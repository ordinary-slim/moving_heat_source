#include "../Problem.h"
#include <Eigen/Core>
#include <Eigen/Sparse>

void Problem::updateForcedDofs() {
  // Has to be last step. Similar treatment to Dirichlet BC
  // Does it have to be last step?
  //
  // ALLOC
  ls->indicesForcedDofs.reserve( ls->indicesForcedDofs.size() + dirichletNodes.size() );
  ls->valuesForcedDofs.reserve( ls->valuesForcedDofs.size() + dirichletNodes.size() );

  vector<int> dirichletIndices = dirichletNodes.getTrueIndices();
  for (int inode : dirichletIndices) {
    ls->indicesForcedDofs.push_back( dofNumbering[ inode ] );
    ls->valuesForcedDofs.push_back( dirichletValues[inode] );
  }
}
