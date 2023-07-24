#include "Mesh.h"
#include "MeshTag.h"
#include "Element.h"
#include <Eigen/Core>
#include <chrono>

namespace mesh
{

Element Mesh::getEntity(int ient, const Connectivity &connectivity, const ReferenceElement *refEl, const ReferenceElement *facetRefEl ) const {

  if (ient < 0) {
    throw std::invalid_argument( "received negative value." );
  }

  Element e;
  e.ient = ient;
  e.setElementType( refEl );
  if (facetRefEl) {
    e.facetRefEl = facetRefEl;
  }

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
  return getEntity( ielem, con_CellPoint, &refCellEl, &refFacetEl );
}
vector<int> mesh::Mesh::findOwnerElements( const Eigen::Vector3d &point ) const {
  vector<int> idxOwnerEl;
  vector<int> potentialOwners;
  //Broad  Phase Search
  //Convert to CGAL point and use bounding boxes tree
  CGAL::Simple_cartesian<double>::Point_3 cgalPoint( point[0], point[1], point[2] );
  tree.all_intersected_primitives( cgalPoint, std::back_inserter( potentialOwners ) );

  //Narrow Phase
  Element cellEl;
  Element facetEl;

  for ( int ielem : potentialOwners ) {
    bool isInside = true;
    cellEl = getElement( ielem );

    std::vector<std::vector<unsigned int>> setsFacetLocalCons = getFacetVertexSets( cellEl.refEl->elementType );

    for ( std::vector<unsigned int> facetLocalCon : setsFacetLocalCons ) {
      facetEl = cellEl.getFacetElement( &facetLocalCon );

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
vector<int> mesh::Mesh::findCollidingElements( const myOBB &obb ) const {
  vector<int> indicesCollidingEls;
  vector<int> potentialCollidingEls;
  //Broad  Phase Search
  //Convert to CGAL aabb and use bounding boxes tree
  auto cgal_aabb = static_cast<inex_K::Iso_cuboid_3>( obb );
  tree.all_intersected_primitives( cgal_aabb, std::back_inserter( potentialCollidingEls ) );

  //Narrow Phase
  for ( int ielem : potentialCollidingEls ) {
    Element cellEl = getElement( ielem );
    if (obb.hasCollided( cellEl )) {
      indicesCollidingEls.push_back( ielem );
    }
  }
  return indicesCollidingEls;
}

void mesh::Mesh::setAABBs() {
  // CGAL AABB tree
  auto begin = std::chrono::steady_clock::now();

  elementAABBs.reserve( nels );

  for (int ielem = 0; ielem < nels; ++ielem) {
    Element e = getElement(ielem);
    elementAABBs.push_back( myAABB( e ) );
  }

  tree.rebuild(elementAABBs.begin(), elementAABBs.end());

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
