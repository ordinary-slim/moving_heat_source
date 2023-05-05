#ifndef SUBMESH
#define SUBMESH
#include <algorithm>
#include "Mesh.h"
#include "MeshTag.h"
#include "MeshTagImpl.h"
#include "Element.h"

namespace mesh {
class Submesh {
  public:
    mesh::Mesh *mesh;//this should be private
                     
    bool hasInactive;
    mesh::MeshTag<int> activeNodes, activeElements, justDeactivatedElements, justActivatedBoundary;
    MeshTag<int> boundaryFacets, boundaryFacetsParentEls;

    Submesh(Mesh *m) :
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
      return mesh->getEntity(ient, connectivity, refEl);
    }
    Element getElement(int ielem) {
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
    void setActivation(const MeshTag<int> &activation);
    void updateActiveNodes();
    void updateActiveElements();
    void updateBeforeActivation();
    void updateAfterActivation();
  private:
    int _dim;
};
}
#endif
