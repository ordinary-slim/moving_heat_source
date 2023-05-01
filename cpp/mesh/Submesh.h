#ifndef SUBMESH
#include <algorithm>
#include "Mesh.h"
#include "MeshTag.h"
#include "Element.h"

namespace mesh {
class Submesh {
  public:
    mesh::Mesh *mesh;//this should be private
                     
    bool hasInactive;
    mesh::MeshTag<int> activeNodes, activeElements;
    mesh::Boundary boundary;

    Submesh() = default;

    Submesh(Mesh *m) {
      mesh = m;
      _dim = mesh->dim;
      activeNodes = mesh::MeshTag<int>(mesh, 0);
      activeElements = mesh::MeshTag<int>(mesh, mesh->dim);
      std::fill(activeNodes.x.begin(), activeNodes.x.end(), 1);
      std::fill(activeElements.x.begin(), activeElements.x.end(), 1);
      boundary = mesh->boundary;
    }

    Boundary findBoundary();
    Element getEntity(int ient, Connectivity &connectivity, ReferenceElement &refEl ) {
      return mesh->getEntity(ient, connectivity, refEl);
    }
    Element getElement(int ielem) {
      return mesh->getElement(ielem);
    }
    Element getBoundaryFacet(int ifacet) {
      // Assumed that ifacet is a boundary facet
      Element e = getEntity( ifacet, mesh->con_FacetPoint, mesh->refFacetEl );
      Element parentEl = getElement( boundary.parentEls[ ifacet ] );
      e.computeNormal( parentEl.getCentroid() );
      return e;
    }

    void setActivation(const MeshTag<int> &activation);
    void updateActiveNodes();
    void updateActiveElements();
    void updateAfterActivation();
  private:
    int _dim;
};
}
#define SUBMESH
#endif
