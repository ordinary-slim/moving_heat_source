#ifndef SUBMESH
#define SUBMESH
#include <algorithm>
#include "Mesh.h"
#include "MeshTag.h"
#include "MeshTagImpl.h"
#include "Element.h"
#include "../linearAlgebra/LinearSystem.h"
#include <Eigen/Sparse>

class Problem;
namespace fem {
  class Function;
}

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
    // For solving projections
    std::shared_ptr<LinearSystem> ls = NULL;
    vector<int> dofNumbering;
    Eigen::SparseMatrix<double, Eigen::RowMajor> *massMat = NULL;// Pointer to lhs of member ls
                     
    bool hasInactive;
    mesh::MeshTag<int> materialTag;// tags each cell with index of its material
    mesh::MeshTag<int> activeNodes, activeElements;
    mesh::MeshTag<int> justDeactivatedElements, justActivatedBoundary;//Are these useful?
    MeshTag<int> boundaryFacets, boundaryFacetsParentEls;

    Domain(Mesh *m, Problem *p);

    void preIterate();
    void assembleMassMatrix();
    void setSpeed(Eigen::Vector3d speedDomain);


    Element getEntity(int ient, Connectivity &connectivity, ReferenceElement &refEl ) {
      return mesh->getEntity(ient, connectivity, &refEl);
    }
    Element getElement(int ielem) const {
      Element e = mesh->getElement(ielem);
      e.imat = materialTag[ ielem ];
      return e;
    }
    Element getBoundaryFacet(int ifacet) {
      // Assumed that ifacet is a boundary facet
      Element e = getEntity( ifacet, mesh->con_FacetPoint, mesh->refFacetEl );
      int iParentElement = boundaryFacetsParentEls[ifacet];
      Element parentEl = getElement(iParentElement);
      e.computeNormal( parentEl.getCentroid() );
      e.imat = materialTag[iParentElement];
      return e;
    }

    int dim() const { return _dim; }

    void computeBoundary();
    int findOwnerElements( const Eigen::Vector3d &point ) const;
    void setMaterialSets(const MeshTag<int> &materialTag);
    void setActivation(const MeshTag<int> &activation);
    void resetActivation();
    void deactivate();
    void updateActiveNodes(const MeshTag<int> *newActiveNodes = NULL);
    void updateActiveElements(const MeshTag<int> *newActiveEls = NULL);
    void updateBeforeActivation();
    void updateAfterActivation();
    void intersect( const MeshTag<int> subdomain);
    void inPlaneRotate( Eigen::Vector3d &center, double angle );
    void invertProjection( Eigen::VectorXd &sol, Eigen::VectorXd &projection );
    fem::Function project( std::function<double(Eigen::Vector3d)> func );//L2 projection onto domain attribute
    MeshTag<int> projectCellTag( const MeshTag<int> &cellTag, const Domain &extDomain );
  private:
    int _dim;
};
  void throwPointOutOfBounds(  const Eigen::Vector3d &point  );
}
#endif
