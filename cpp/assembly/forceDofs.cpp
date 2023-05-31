#include "../Problem.h"
#include <Eigen/Core>
#include <Eigen/Sparse>

void Problem::forceDofs() {
  // ALLOC
  ls->indicesForcedDofs.reserve( ls->indicesForcedDofs.size() + forcedDofs.size() );
  ls->valuesForcedDofs.reserve( ls->valuesForcedDofs.size() + forcedDofs.size() );

  vector<int> forcedDofsIndices = forcedDofs.getTrueIndices();
  for (int inode : forcedDofsIndices) {
    ls->indicesForcedDofs.push_back( dofNumbering[ inode ] );
    ls->valuesForcedDofs.push_back( forcedDofsValues[inode] );
  }
}
