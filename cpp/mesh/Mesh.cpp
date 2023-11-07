#include "Mesh.h"
#include "Element.h"
#include <Eigen/Core>
#include <cmath>
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
vector<int> Mesh::findOwnerElements( const Eigen::Vector3d &point ) const {
  vector<int> indicesOwnerElements;
  vector<int> potentialOwners;
  //Broad  Phase Search
  //Convert to CGAL point and use bounding boxes tree
  CGAL::Simple_cartesian<double>::Point_3 cgalPoint( point[0], point[1], point[2] );
  looseTree.all_intersected_primitives( cgalPoint, std::back_inserter( potentialOwners ) );

  //Narrow Phase
  Element cellEl;
  Element facetEl;

  indicesOwnerElements.reserve( potentialOwners.size() );
  for ( int ielem : potentialOwners ) {
    bool isInside = true;
    cellEl = getElement( ielem );

    std::vector<std::vector<unsigned int>> setsFacetLocalCons = getFacetVertexSets( cellEl.refEl->elementType );

    for ( std::vector<unsigned int>& facetLocalCon : setsFacetLocalCons ) {
      facetEl = cellEl.getFacetElement( &facetLocalCon );

      double projection =  facetEl.normal.dot( point - facetEl.centroid );
      if ( projection > +toleranceSearches ) {
        isInside = false;
        break;
      }
    }
    if (isInside) {
      indicesOwnerElements.push_back( ielem );
    }
  }
  return indicesOwnerElements;
}
vector<int> Mesh::findCollidingElements( const MyOBB &obb ) const {
  vector<int> indicesCollidingEls;
  vector<int> potentialCollidingEls;
  //Broad  Phase Search
  //Convert to CGAL aabb and use bounding boxes tree
  auto cgal_aabb = static_cast<inex_K::Iso_cuboid_3>( obb );
  looseTree.all_intersected_primitives( cgal_aabb, std::back_inserter( potentialCollidingEls ) );

  //Narrow Phase
  indicesCollidingEls.reserve( potentialCollidingEls.size() );
  for ( int ielem : potentialCollidingEls ) {
    Element cellEl = getElement( ielem );
    if (obb.hasCollided( cellEl )) {
      indicesCollidingEls.push_back( ielem );
    }
  }
  return indicesCollidingEls;
}

vector<int> Mesh::findCollidingElements( const MyAABB &aabb ) const {
  vector<int> indicesCollidingEls;
  //Single search against tight AABBs
  auto cgal_aabb = static_cast<inex_K::Iso_cuboid_3>( aabb );
  tightTree.all_intersected_primitives( cgal_aabb, std::back_inserter( indicesCollidingEls ) );
  return indicesCollidingEls;
}


vector<int> Mesh::findCollidingElements( const Eigen::Vector3d &center, const double R) const {
  vector<int> indicesCollidingEls;
  vector<int> potentialCollidingEls;
  //Broad  Phase Search
  //Convert to CGAL aabb and use bounding boxes tree
  double minX = center[0] - R, maxX = center[0] + R,
         minY = center[1] - R, maxY = center[1] + R,
         minZ = center[2] - R, maxZ = center[2] + R;
  auto cgal_aabb = inex_K::Iso_cuboid_3( minX, minY, minZ, maxX, maxY, maxZ );
  looseTree.all_intersected_primitives( cgal_aabb, std::back_inserter( potentialCollidingEls ) );

  // Narrow phase
  indicesCollidingEls.reserve( potentialCollidingEls.size() );
  for ( int ielem : potentialCollidingEls ) {
    Element cellEl = getElement( ielem );
    // Compute distance to center of ball
    double distance = (cellEl.centroid - center).norm();
    // Compare to cutoff
    if ( distance <= R ) {
      indicesCollidingEls.push_back( ielem );
    }
  }
  return indicesCollidingEls;
}

void Mesh::buildAABBTrees() {
  // CGAL AABB tree
  looseElementAABBs.resize( nels );
  tightElementAABBs.resize( nels );
  updateAABBTrees();
}

void Mesh::updateAABBTrees() {
  for (int ielem = 0; ielem < nels; ++ielem) {
    Element e = getElement(ielem);
    looseElementAABBs[ielem] = MyAABB( e, toleranceSearches, 0.05 );
    tightElementAABBs[ielem] = MyAABB( e, 0.0, 0.0 );
  }
  looseTree.rebuild(looseElementAABBs.begin(), looseElementAABBs.end());
  tightTree.rebuild(tightElementAABBs.begin(), tightElementAABBs.end());
}


void inPlaneRotate( Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> &points, 
    const Eigen::Vector3d &center, double angle ) {
  /*
   * points: Array of points to be rotated, row by row
   * Center: Point along axis of rotation
   * angle: In radians
   */

  // Build rotation matrix
  double c = std::cos(angle);
  double s = std::sin(angle);

  Eigen::Matrix3d R {
    {+c, -s,  0},
    {+s, +c,  0},
    { 0,  0,  1}
  };
  // Loop over points
  for (int ipoin = 0; ipoin < points.rows(); ++ipoin) {
    // Rotate each point
    points.row( ipoin ) = center + R * ( points.row( ipoin ).transpose() - center );
  }
  // Round to N decimal places
  int decimal_places = 7;
  const double multiplier = std::pow(10.0, decimal_places);
  points.unaryExpr( [multiplier]( double d ) { return std::round( d * multiplier ) / multiplier; } );
}

}
