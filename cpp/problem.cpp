#include "problem.h"
#include "FEMFunction.h"
#include "mesh/Mesh.h"
#include <Eigen/Core>

void Problem::updateFRFpos() {
  // Pre-iteration operations

  // update positions in no advection RF
  // done in pre iterate because activation is also done here
  shiftFRF += -dt * advectionSpeed;
  for (int inode=0; inode < mesh.nnodes; inode++){
    mesh.posFRF.row( inode ) += -dt * advectionSpeed;
  }
}

void Problem::getFromExternal(mesh::Mesh &extMesh, FEMFunction &extFEMFunc, Eigen::Vector3d shiftExtFRF ){
  Eigen::Vector3d posExt;
  std::fill(unknown.values.begin(), unknown.values.end(), 0.0);
  for (int inode = 0; inode < mesh.nnodes; inode++) {
    posExt = mesh.pos.row(inode) + (shiftFRF - shiftExtFRF).transpose();
    unknown.values[inode] = extFEMFunc.interpolate( posExt );
  }
}
