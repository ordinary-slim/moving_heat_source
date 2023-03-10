#ifndef ELEMENT
#include <iostream>
#include <vector>
#include <Eigen/Core>
#include <algorithm>
#include "RefElement.h"

namespace mesh
{
class Element {
  public:
    int nnodes, ngpoints;
    Eigen::MatrixX3d pos, gpos;
    Eigen::VectorXi  con;
    std::vector<double> gpweight;
    double vol;
    Eigen::Vector3d centroid;
    Eigen::Vector3d normal;
    int dim;
    ElementType elementType;
    std::vector<std::vector<double>> BaseGpVals;
    std::vector<std::vector<Eigen::Vector3d>> GradBaseGpVals;
    bool openIntegration = false;//default closed integration

    ReferenceElement      *refEl;
    // x : loc coordinate; xi : reference coordinate
    // x = ref2locShift + ref2locMatrix · xi
    // xi = loc2refShift + loc2refMatrix · x
    Eigen::Vector3d ref2locShift, loc2refShift;
    Eigen::Matrix3d ref2locMatrix, loc2refMatrix;

    void computeNodalValues_Base(){
      //COMMON
      //TODO: Take all of this from reference element
      BaseGpVals = refEl->BaseGpVals;
      gpweight = refEl->gpweight;
      openIntegration = refEl->openIntegration;

      if (openIntegration) {
        for (int igpoin = 0; igpoin < ngpoints; ++igpoin) {
          gpos.row(igpoin) = map_ref2loc( refEl->gpos.row(igpoin) );
        }
      } else {
        gpos = pos;
      }
    }

    void computeNodalValues_GradBase() {
      computeDerivatives();
    }

    void allocate() {
      // ALLOCATIONS
      pos.resize( nnodes, 3 );
      gpos.resize( ngpoints, 3 );
      // BaseFun
      BaseGpVals.resize( nnodes );
      for (int inode = 0; inode < nnodes; ++inode) {
        BaseGpVals[inode].resize( ngpoints );
      }
      // Quadrature weights
      gpweight.resize( ngpoints );
      // GradBaseFun
      GradBaseGpVals.resize( nnodes );
      for (int inode = 0; inode < nnodes; inode++) {
        GradBaseGpVals[inode].resize( ngpoints );
        for (int jgp = 0; jgp < ngpoints; jgp++) {
          GradBaseGpVals[inode][jgp].setZero();
        }
      }
    }

    void setElementType( ReferenceElement &target_refEl ) {
      refEl = &target_refEl;
      elementType = refEl->elementType;
      nnodes = refEl->nnodes;
      ngpoints = refEl->ngpoints;
      dim = refEl->dim;
    }

    void computeLocRefMappings() {
      // Build helper matrices
      Eigen::Matrix3d X;
      // Compute
      X.setZero();
      for (int inode = 0; inode<dim; inode++) {
        X.row( inode )  = pos.row( inode+1 ) - pos.row( 0 );
      }
      X.transposeInPlace();
      // Fill missing dims
      for (int idim = dim; idim < 3; idim++) {
        X(idim, idim)  = 1.0;
      }
      
      // Compute local <-> reference mappings
      ref2locMatrix = X * refEl->XI_inverse;
      ref2locShift  = pos.row(0).transpose() - ref2locMatrix * refEl->pos.row(0).transpose() ;
      loc2refMatrix = ref2locMatrix.inverse();
      loc2refShift  = refEl->pos.row(0).transpose() - loc2refMatrix * pos.row(0).transpose();
    }

    Eigen::Vector3d map_ref2loc( Eigen::Vector3d xi ) {
      return ref2locShift + ref2locMatrix * xi;
    }

    Eigen::Vector3d map_loc2ref( Eigen::Vector3d x ) {
      return loc2refShift + loc2refMatrix * x;
    }

    void computeDerivatives() {
      Eigen::Matrix3d loc2refMatrix_T = loc2refMatrix.transpose();

      // Compute volume
      vol = refEl->vol * ref2locMatrix.determinant();

      for (int inode = 0; inode < nnodes; inode++) {
        for (int igp = 0; igp < ngpoints; igp++) {
          GradBaseGpVals[inode][igp] = loc2refMatrix_T*(refEl->GradBaseGpVals[inode][igp]);
        }
      }
    }


    void computeCentroid() {
      centroid.setZero();
      for (int inode = 0; inode < nnodes; ++inode){
        centroid += pos.row(inode);
      }
      centroid /= nnodes;
    }

    Eigen::Vector3d getCentroid() {
      return centroid;
    }

    void computeNormal( Eigen::Vector3d parentCentroid ) {
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

    Element getFacetElement( Eigen::VectorXi vertices, ReferenceElement &facetRefEl ) {
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
      e.computeNormal( getCentroid() );
      return e;
    }

    Eigen::VectorXd evaluateShaFuns( Eigen::Vector3d pos ) {
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

    double getSizeAlongVector( Eigen::Vector3d vector );
};
}
#define ELEMENT
#endif
