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
      activeNodes(mesh::MeshTag<int>(m, 0)),
      activeElements(mesh::MeshTag<int>(m, m->dim)),
      justDeactivatedElements(mesh::MeshTag<int>(m, m->dim)),
      justActivatedBoundary(mesh::MeshTag<int>(m, m->dim-1)),
      boundaryFacets(mesh::MeshTag<int>(m, m->dim-1)),
      boundaryFacetsParentEls(mesh::MeshTag<int>(m, m->dim-1))
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
      int idxParentEl = mesh->con_FacetCell.con[ifacet][0];
      Element parentEl = getElement( idxParentEl );
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
