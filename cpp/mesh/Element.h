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
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>
      pos, gpos;
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

    Eigen::Vector3d getCentroid() {
      return centroid;
    }

    void computeCentroid();
    void computeNormal( Eigen::Vector3d parentCentroid );
    Element getFacetElement( Eigen::VectorXi vertices, ReferenceElement &facetRefEl );
    Eigen::VectorXd evaluateShaFuns( Eigen::Vector3d pos );
    double getSizeAlongVector( Eigen::Vector3d vector );
};
}
#define ELEMENT
#endif
