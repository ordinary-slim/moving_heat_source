#include "mesh.h"
#include "element.h"
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
}
