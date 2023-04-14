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
    e.pos.row(inode) = pos.row(e.con(inode));
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
  Eigen::VectorXi facets;
  Element cellEl;
  Element facetEl;

  for ( int ielem : potentialOwners ) {
    bool isInside = true;
    cellEl = getElement( ielem );
    facets = con_CellFacet.getLocalCon( ielem );
    for ( int ifacet : facets ) {
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

void mesh::Mesh::setActiveElements(vector<int> inputActiveElements ) {
  activeElements = inputActiveElements;
  //Update activeNodes
  fill( activeNodes.begin(), activeNodes.end(), 0 );
  for (int ielem = 0; ielem < nels; ielem++) {
    if (activeElements[ielem] == 1) {
      Eigen::VectorXi locCon = con_CellPoint.getLocalCon( ielem );
      //set to 1 nodes who belong to element
      for (int locInode = 0; locInode < locCon.size(); locInode++){
        activeNodes[ locCon[locInode] ] = 1;
      }
    }
  }
  if (std::find( activeElements.begin(), activeElements.end(), 0)
      != activeElements.end() ) {
    hasInactive = true;
    findBoundary();
  } else {
    hasInactive = false;
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
    for (int iel = 0; iel < con_FacetCell.con.cols(); ++iel) {
      if (con_FacetCell.con(ifacet, iel) == -1) { break; }
      
      if (activeElements[ con_FacetCell.con(ifacet, iel) ]) {
        lastVisitedActiveEl = con_FacetCell.con(ifacet, iel);
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
