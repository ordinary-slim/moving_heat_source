#ifndef SUBMESH
#define SUBMESH
#include <algorithm>
#include "Mesh.h"
#include "MeshTag.h"
#include "MeshTagImpl.h"
#include "Element.h"
#include <Eigen/Sparse>

namespace mesh {
/*
 * Wrapper around Mesh object
 * Defines a subset of the Mesh
 */
class ActiveMesh {
  public:
    mesh::Mesh *mesh;//this should be private

    // Mass matrix
    Eigen::SparseMatrix<double> massMat;
    vector<Eigen::Triplet<double>> massCoeffs;
                     
    bool hasInactive;
    mesh::MeshTag<int> activeNodes, activeElements;
    mesh::MeshTag<int> justDeactivatedElements, justActivatedBoundary;//Are these useful?
    MeshTag<int> boundaryFacets, boundaryFacetsParentEls;

    ActiveMesh(Mesh *m) :
      activeNodes(mesh::MeshTag<int>(m, 0, 1)),
      activeElements(mesh::MeshTag<int>(m, m->dim, 1)),
      justDeactivatedElements(mesh::MeshTag<int>(m, m->dim, 0)),
      justActivatedBoundary(mesh::MeshTag<int>(m, m->dim-1, 0)),
      boundaryFacets(mesh::MeshTag<int>(m, m->dim-1)),
      boundaryFacetsParentEls(mesh::MeshTag<int>(m, m->dim-1))
    {
      mesh = m;
      _dim = mesh->dim;
      computeBoundary();
    }

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
  private:
    int _dim;
};
}
#endif
