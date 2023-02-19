#include "problem.h"
#include "FEMFunction.h"
#include "mesh/Mesh.h"

void Problem::getFromExternal(mesh::Mesh &extMesh, FEMFunction &extFEMFunc){
  std::fill(unknown.values.begin(), unknown.values.end(), 0.0);
  for (int inode = 0; inode < mesh.nnodes; inode++) {
    unknown.values[inode] = extFEMFunc.interpolate( mesh.pos.row(inode) );
  }
}
