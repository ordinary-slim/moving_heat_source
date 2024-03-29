#ifndef ELEMENT
#include <iostream>
#include <vector>
#include <Eigen/Core>
#include <algorithm>
#include "RefElement.h"

typedef Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> Dense3ColMat;

namespace mesh
{
class Element {
  public:
    int nnodes, ngpoints;
    Dense3ColMat pos, gpos;
    const std::vector<int>*  con;
    std::vector<double> gpweight;
    double vol;
    Eigen::Vector3d centroid;
    Eigen::Vector3d normal;
    int dim;
    ElementType elementType;
    std::vector<std::vector<double>> BaseGpVals;
    std::vector<std::vector<Eigen::Vector3d>> GradBaseGpVals;
    bool openIntegration = false;//default closed integration

    const ReferenceElement      *refEl;
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

    void setElementType( const ReferenceElement &target_refEl ) {
      refEl = &target_refEl;
      elementType = refEl->elementType;
      nnodes = refEl->nnodes;
      ngpoints = refEl->ngpoints;
      dim = refEl->dim;
    }

    void computeLocRefMappings();

    void computeDerivatives();

    Eigen::Vector3d map_ref2loc( Eigen::Vector3d xi ) {
      return ref2locShift + ref2locMatrix * xi;
    }

    Eigen::Vector3d map_loc2ref( Eigen::Vector3d x ) {
      return loc2refShift + loc2refMatrix * x;
    }

    Eigen::Vector3d getCentroid() const {
      return centroid;
    }

    void computeCentroid();
    void computeNormal( Eigen::Vector3d parentCentroid );
    Element getFacetElement( const std::vector<int>* vertices, ReferenceElement &facetRefEl ) const;
    Eigen::VectorXd evalShaFuns( Eigen::Vector3d pos );
    Dense3ColMat evaluateGradShaFuns( Eigen::Vector3d pos );
    double getSizeAlongVector( Eigen::Vector3d vector );
};
}
#define ELEMENT
#endif
