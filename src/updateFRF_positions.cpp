#include "problem.h"

void Problem::updateFRF_positions() {
  // Pre-iteration operations

  // update positions in no advection RF
  // done in pre iterate because activation is also done here
  for (int inode=0; inode < mesh.nnodes; inode++){
    mesh.pos_noAdv.row( inode ) += -dt * advectionSpeed;
  }
}
