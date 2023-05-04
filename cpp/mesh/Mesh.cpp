#include "Mesh.h"
#include "Element.h"
#include <Eigen/Core>
#include <chrono>

namespace mesh
{

Element Mesh::getEntity(int ient, Connectivity &connectivity, ReferenceElement &refEl ) {

  if (ient < 0) {
    throw std::invalid_argument( "received negative value." );
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
  Element parentEl = getElement( boundary.parentEls[ ifacet ] );
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
  const vector<int>* facets;
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

Boundary mesh::Mesh::findBoundary() {
  /*
   * Build array of indices of boundary facets
  */
  Boundary b;
  b.facets.clear();
  b.parentEls.clear();
  b.parentEls.resize( con_FacetCell.nels_oDim );
  std::fill( b.parentEls.begin(), b.parentEls.end(), -1 );
  for (int ifacet = 0; ifacet < con_FacetCell.nels_oDim; ++ifacet) {
    const vector<int>* incidentElements = con_FacetCell.getLocalCon( ifacet );
    if (incidentElements->size()==1) {
      b.facets.push_back( ifacet );
      b.parentEls[ifacet] = (*incidentElements)[0];
    }
  }
  return b;
}
}
