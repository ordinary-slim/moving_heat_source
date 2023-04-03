#include "Element.h"
#include <algorithm>

namespace mesh
{
void mesh::Element::computeLocRefMappings() {
  // Build helper matrices
  Eigen::Matrix3d X;
  // Compute
  X.setZero();
  for (int inode = 0; inode<dim; inode++) {
    X.row( inode )  = pos.row( inode+1 ) - pos.row( 0 );
  }

  // Obtain second vector
  if (dim < 2) {
    double tol = 1e-7;
    //Eigen::Vector3d e1 = X.row( 0 );
    Eigen::Vector3d e2;
    // Try cross product with canonical basis
    // Looking for a non-zero vector
    e2 << 0.0, X(0,2), -X(0,1);// x1 x e1
    if (e2.norm() < tol) {
      e2 << -X(0,2), 0.0, X(0,0);// x1 x e2
      if (e2.norm() < tol) {
        e2 << X(0,1), -X(0,0), 0.0;// x1 x e3
      }
    }
    e2 /= e2.norm();

    X.row( 1 ) = e2;
  }
  if (dim < 3) {
    X.row( 2 ) = X.row(0).cross( X.row(1) );// e1 x e2 = e3 <=> right handed
    X.row( 2 ) /= X.row(2).norm();
  }
  X.transposeInPlace();
  
  // Compute local <-> reference mappings
  ref2locMatrix = X * refEl->XI_inverse;
  ref2locShift  = pos.row(0).transpose() - ref2locMatrix * refEl->pos.row(0).transpose() ;
  loc2refMatrix = ref2locMatrix.inverse();
  loc2refShift  = refEl->pos.row(0).transpose() - loc2refMatrix * pos.row(0).transpose();
}

void mesh::Element::computeDerivatives() {
  Eigen::Matrix3d loc2refMatrix_T = loc2refMatrix.transpose();

  // Compute volume
  vol = refEl->vol * ref2locMatrix.determinant();

  for (int inode = 0; inode < nnodes; inode++) {
    for (int igp = 0; igp < ngpoints; igp++) {
      GradBaseGpVals[inode][igp] = loc2refMatrix_T*(refEl->GradBaseGpVals[inode][igp]);
    }
  }
}

Element mesh::Element::getFacetElement( Eigen::VectorXi vertices, ReferenceElement &facetRefEl ) {
  Element e;
  e.setElementType( facetRefEl );
  e.allocate();

  // positions
  int locInode;
  int counter = 0;
  for (int globInode : vertices ) {
    locInode = std::find(con.begin(), con.end(), globInode) - con.begin();
    e.pos.row( counter ) = pos.row( locInode );
    counter++;
  }
  // computations
  e.computeCentroid();
  e.computeNormal( getCentroid() );
  return e;
}

Eigen::VectorXd mesh::Element::evaluateShaFuns( Eigen::Vector3d pos ) {
  // Evaluate shape funcs at a point
  Eigen::VectorXd shaFunsVals( nnodes );
  
  // Map to reference element
  Eigen::Vector3d xi = map_loc2ref( pos );

  // Evaluate shafuns in reference element
  for (int inode = 0; inode < nnodes; inode++) {
    shaFunsVals(inode) = refEl->shapeFuns[inode]( xi );
  }

  return shaFunsVals;
}

void mesh::Element::computeNormal( Eigen::Vector3d parentCentroid ) {
  normal.setZero();
  switch (dim) {
    case 0: {//point1
      normal(0) = 1;
      break;
    }
    case 1: {//line2
      Eigen::Vector3d lineVec = pos.row(1) - pos.row(0);
      normal(0) = -lineVec( 1 );
      normal(1) = +lineVec( 0 );
      normal /= lineVec.norm();
      break;
    }
    case 2: {
      printf("Normal computation for 2D els not implemented yet\n");
      exit(EXIT_FAILURE);
    }
    default: {
      printf("This normal computation is not implemented yet.");
      exit(EXIT_FAILURE);
    }
  }

  // Test for sign
  Eigen::Vector3d testVector = parentCentroid - centroid;
  if ( testVector.dot( normal ) > 0 ) normal = -normal;
}

void mesh::Element::computeCentroid() {
  centroid.setZero();
  for (int inode = 0; inode < nnodes; ++inode){
    centroid += pos.row(inode);
  }
  centroid /= nnodes;
}

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
