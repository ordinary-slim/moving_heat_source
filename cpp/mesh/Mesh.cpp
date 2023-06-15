#include "Mesh.h"
#include "MeshTag.h"
#include "Element.h"
#include <Eigen/Core>
#include <chrono>

namespace mesh
{

Element Mesh::getEntity(int ient, const Connectivity &connectivity, const ReferenceElement &refEl ) const {

  if (ient < 0) {
    throw std::invalid_argument( "received negative value." );
  }

  Element e;
  e.ient = ient;
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

Element Mesh::getElement(int ielem) const {
  return getEntity( ielem, con_CellPoint, refCellEl );
}
vector<int> mesh::Mesh::findOwnerElement( const Eigen::Vector3d &point ) {
  vector<int> idxOwnerEl;
  vector<int> potentialOwners;
  //Broad  Phase Search
  for (int ielem = 0; ielem < nels; ++ielem) {
    if ( elementAABBs[ielem].isPointInside( point ) ) {
      potentialOwners.push_back( ielem );
    }
  }

  //Narrow Phase
  const vector<unsigned int>* facets;
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
      idxOwnerEl.push_back( ielem );
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
MeshTag<int> mark( const Mesh &mesh, int dim, const vector<int> &indices ) {
  vector<int> values = vector<int>( mesh.getNumEntities( dim ), 0 );
  for (int index : indices) {
    values[index] = 1;
  }
  return MeshTag<int>( &mesh, dim, values );
}
}
