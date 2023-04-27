#include "Mesh.h"
#include "Element.h"
#include <Eigen/Core>
#include <chrono>

namespace mesh
{

Element Mesh::getEntity(int ient, Connectivity &connectivity, ReferenceElement &refEl ) {

  if (ient < 0) {
    cout << "Bad Mesh::getEntity!" << endl;
    exit(-1);
  }

  Element e;
  e.setElementType( refEl );

  e.allocate();

  // set connectivity
  e.con = connectivity.getLocalCon( ient );
  // set pos
  for (int inode=0; inode < e.nnodes; inode++) {
    e.pos.row(inode) = pos.row((*e.con)[inode]);
  }

  // COMPUTATIONS
  e.computeCentroid();
  e.computeLocRefMappings();
  e.computeNodalValues_Base();//COMMON BETWEEN ELS
  e.computeNodalValues_GradBase();//UNCOMMON
  return e;
}

Element Mesh::getElement(int ielem) {
  return getEntity( ielem, con_CellPoint, refCellEl );
}
Element Mesh::getBoundaryFacet(int ifacet) {
  // Assumed that ifacet is a boundary facet
  Element e = getEntity( ifacet, con_FacetPoint, refFacetEl );
  Element parentEl = getElement( boundaryFacetsParentEl[ ifacet ] );
  e.computeNormal( parentEl.getCentroid() );
  return e;
}

int mesh::Mesh::findOwnerElement( Eigen::Vector3d point ) {
  int idxOwnerEl = -1;
  vector<int> potentialOwners;
  //Broad  Phase
  for (int ielem = 0; ielem < nels; ++ielem) {
    if ( elementAABBs[ielem].isPointInside( point ) ) {
      potentialOwners.push_back( ielem );
    }
  }

  //Narrow Phase
  vector<int>* facets;
  Element cellEl;
  Element facetEl;

  for ( int ielem : potentialOwners ) {
    bool isInside = true;
    cellEl = getElement( ielem );
    facets = con_CellFacet.getLocalCon( ielem );
    for ( int ifacet : *facets ) {
      facetEl = cellEl.getFacetElement( con_FacetPoint.getLocalCon(ifacet),
                                        refFacetEl);

      if ( facetEl.normal.dot( point - facetEl.centroid )  > 0 ) {
        isInside = false;
        break;
      }
    }
    if (isInside) {
      idxOwnerEl = ielem;
      break;
    }
  }
  return idxOwnerEl;
}
void mesh::Mesh::setAABBs() {

  auto begin = std::chrono::steady_clock::now();

  elementAABBs.resize( nels );

  Element e;

  for (int ielem = 0; ielem < nels; ++ielem) {
    e = getElement(ielem);
    elementAABBs[ielem] = AABB( e );
  }

  auto end = std::chrono::steady_clock::now();
  std::cout << "Building AABBs took " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
}

void mesh::Mesh::updateActiveNodes() {
  /*
   * Update activeNodes after a change in activeElements
   * If a node belongs to an active element, set it to active.
   */
  vector<int>* incidentElements;
  for (int inode = 0; inode < nnodes; ++inode) {
    activeNodes[inode] = 0;
    incidentElements = con_PointCell.getLocalCon( inode );
    for ( auto p_ielem = incidentElements->begin(); (p_ielem != incidentElements->end())&&(*p_ielem != -1); ++p_ielem ) {
      if (activeElements[*p_ielem] == 1) {
        activeNodes[inode] = 1;
        break;
      }
    }
  }
}

void mesh::Mesh::updateActiveElements() {
  /*
   * Update activeElements after a change in activeNodes
   * If all the nodes of an element are active, activate it.
   */
  vector<int>* incidentNodes;
  bool allNodesActive;
  for (int ielem = 0; ielem < nels; ++ielem) {
    activeElements[ielem] = 0;
    incidentNodes = con_CellPoint.getLocalCon( ielem );
    allNodesActive = true;
    for ( auto p_inode = incidentNodes->begin(); (p_inode != incidentNodes->end())&&(*p_inode != -1); ++p_inode ) {

      if (activeNodes[*p_inode] == 0) {
        allNodesActive = false;
        break;
      }
    }
    if (allNodesActive) {
      activeElements[ielem] = 1;
    }
  }
}

bool checkHasInactive( const vector<int> &activeElements ) {
  return (std::find( activeElements.begin(), activeElements.end(), 0) != activeElements.end() );
}

void mesh::Mesh::setActiveElements(const vector<int> &otherActiveElements ) {
  activeElements = otherActiveElements;
  updateActiveNodes();
  hasInactive = checkHasInactive(activeElements);
  if (hasInactive) {
    findBoundary();
  }
}

void mesh::Mesh::setActiveNodes(const vector<int> &otherActiveNodes ) {
  activeNodes = otherActiveNodes;
  updateActiveElements();
  hasInactive = checkHasInactive(activeElements);
  if (hasInactive) {
    findBoundary();
  }
}

void mesh::Mesh::findBoundary() {
  /*
   * Build array of indices of boundary facets
   * Check second element of con and decide
  */
  boundaryFacets.clear();
  boundaryFacetsParentEl.clear();
  boundaryFacetsParentEl.resize( con_FacetCell.nels_oDim );
  std::fill( boundaryFacetsParentEl.begin(), boundaryFacetsParentEl.end(), -1 );
  int activeElsPerFacet;
  int lastVisitedActiveEl;
  for (int ifacet = 0; ifacet < con_FacetCell.nels_oDim; ++ifacet) {
    activeElsPerFacet = 0;
    vector<int>* incidentElements = con_FacetCell.getLocalCon( ifacet );
    for ( auto p_ielem = incidentElements->begin(); (p_ielem != incidentElements->end())&&(*p_ielem != -1); ++p_ielem ) {
      
      if (activeElements[ *p_ielem ]) {
        lastVisitedActiveEl = *p_ielem;
        ++activeElsPerFacet;
      }
    }
    if (activeElsPerFacet==1) {
      boundaryFacets.push_back( ifacet );
      boundaryFacetsParentEl[ifacet] = lastVisitedActiveEl;
    }
  }
  return;
}

}
