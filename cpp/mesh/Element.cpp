#include "Element.h"
#include <algorithm>

namespace mesh
{
double mesh::Element::getSizeAlongVector( Eigen::Vector3d vector ){
  // Relies on passing vector by copy

  vector /= vector.norm();

  double* projections = new double[nnodes];

  // Project each point along vector
  for (int inode = 0; inode < nnodes; ++inode) {
    projections[inode] = vector.dot(pos.row( inode ));
  }
  // Compute max and min, substract, return
  double maxProjec = projections[0];
  double minProjec = projections[0];
  for (int inode = 1; inode < nnodes; ++inode) {
    maxProjec = std::max( maxProjec, projections[inode] );
    minProjec = std::min( minProjec, projections[inode] );
  }


  delete[] projections;
  return (maxProjec-minProjec);
}
}
