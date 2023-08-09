#ifndef SUBMESH
#define SUBMESH
#include <algorithm>
#include "Mesh.h"
#include "MeshTag.h"
#include "MeshTagImpl.h"
#include "Element.h"
#include <Eigen/Sparse>

class Problem;

namespace mesh {
/*
 * Wrapper around Mesh object
 * Defines a subset of the Mesh
 * and keeps track of mesh motion
 * to communicate with problem in Fixed /Laboratory Reference Frame
 */
class Domain {
  private:
    const Problem *problem;//observer
  public:
    mesh::Mesh *mesh;//this should be private
                     //
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>
      posLab; // node positions in fixed / lab reference frame
              // position of point with idx i = row(i)
    Eigen::Vector3d translationLab = Eigen::Vector3d::Zero();//pos + translation = posFRF
    Eigen::Vector3d speedDomain = Eigen::Vector3d::Zero();//domain speed with
                                      //respect to Fixed / Lab reference frame

    // Mass matrix
    Eigen::SparseMatrix<double> massMat;
    vector<Eigen::Triplet<double>> massCoeffs;
                     
    bool hasInactive;
    mesh::MeshTag<int> activeNodes, activeElements;
    mesh::MeshTag<int> justDeactivatedElements, justActivatedBoundary;//Are these useful?
    MeshTag<int> boundaryFacets, boundaryFacetsParentEls;

    Domain(Mesh *m, Problem *p);

    void preIterate();
    void setSpeed(Eigen::Vector3d speedDomain);


    Element getEntity(int ient, Connectivity &connectivity, ReferenceElement &refEl ) {
      return mesh->getEntity(ient, connectivity, &refEl);
    }
    Element getElement(int ielem) const {
      return mesh->getElement(ielem);
    }
    Element getBoundaryFacet(int ifacet) {
      // Assumed that ifacet is a boundary facet
      Element e = getEntity( ifacet, mesh->con_FacetPoint, mesh->refFacetEl );
      Element parentEl = getElement( boundaryFacetsParentEls[ifacet] );
      e.computeNormal( parentEl.getCentroid() );
      return e;
    }

    int dim() { return _dim; }

    void computeBoundary();
    int findOwnerElements( const Eigen::Vector3d &point ) const;
    void setActivation(const MeshTag<int> &activation);
    void resetActivation();
    void deactivate();
    void updateActiveNodes(const MeshTag<int> *newActiveNodes = NULL);
    void updateActiveElements(const MeshTag<int> *newActiveEls = NULL);
    void updateBeforeActivation();
    void updateAfterActivation();
    void intersect( const MeshTag<int> subdomain);
    void inPlaneRotate( Eigen::Vector3d &center, double angle );
  private:
    int _dim;
};
}
#endif
