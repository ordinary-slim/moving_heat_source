#include "Mesh.h"
#include "Element.h"
#include <Eigen/Core>
#include <chrono>

namespace mesh
{
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
  } else {
    hasInactive = false;
  }
}

}
