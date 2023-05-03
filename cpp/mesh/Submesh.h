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
    mesh::MeshTag<int> activeNodes, activeElements, justDeactivatedElements, justActivatedBoundary;
    mesh::Boundary boundary;

    Submesh(Mesh *m) :
      activeNodes(mesh::MeshTag<int>(m, 0)),
      activeElements(mesh::MeshTag<int>(m, m->dim)),
      justDeactivatedElements(mesh::MeshTag<int>(m, m->dim)),
      justActivatedBoundary(mesh::MeshTag<int>(m, m->dim-1))
    {
      mesh = m;
      _dim = mesh->dim;
      activeElements = mesh::MeshTag<int>(mesh, mesh->dim);
      justDeactivatedElements = mesh::MeshTag<int>(mesh, mesh->dim);
      justActivatedBoundary = mesh::MeshTag<int>(mesh, mesh->dim-1);
      std::fill(activeNodes.x.begin(), activeNodes.x.end(), 1);
      std::fill(activeElements.x.begin(), activeElements.x.end(), 1);
      std::fill(justDeactivatedElements.x.begin(), justDeactivatedElements.x.end(), 0);
      std::fill(justActivatedBoundary.x.begin(), justActivatedBoundary.x.end(), 0);
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

    int dim() { return _dim; }

    void setActivation(const MeshTag<int> &activation);
    void updateActiveNodes();
    void updateActiveElements();
    void updateBeforeActivation();
    void updateAfterActivation();
  private:
    int _dim;
};
}
#define SUBMESH
#endif
