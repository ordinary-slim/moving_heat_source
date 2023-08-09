#ifndef MESH
#define MESH
#include <iostream>
#include <vector>
#include <Eigen/Core>
#include "Connectivity.h"
#include "Element.h"
#include "RefElement.h"
#include "ElementTypes.h"
#include "../../external/pybind11/include/pybind11/pybind11.h"
#include "MeshTag.h"
#include "cgal_interface.h"

using namespace std;
using AABB_tree = CGAL::AABB_tree<CGAL::AABB_traits<CGAL::Simple_cartesian<double>,
      myBboxPrimitive>>;

namespace mesh
{
class Mesh {
  public:
    // Copy constructor
    //Mesh(const Mesh& mesh) = default;
    // Construct from python dict with points, cells and cell_type
    Mesh(const pybind11::dict &input);
    int dim;
    int nels, nnodes, nnodes_per_el;
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>
      pos; // nodal positions of point with idx i = row(i)
    Connectivity  con_CellPoint;
    Connectivity  con_PointCell;
    Connectivity  con_CellCell;
    Connectivity  con_FacetPoint;
    Connectivity  con_FacetCell;
    Connectivity  con_CellFacet;
    vector<ElementType> elementTypes;
    ReferenceElement refCellEl;
    ReferenceElement refFacetEl;
    // Fast spatial search
    vector<myAABB> elementAABBs;
    AABB_tree tree;

    Element getEntity(int ient, const Connectivity &connectivity, const ReferenceElement *refEl, const ReferenceElement *facetRefEl = NULL ) const;
    Element getElement(int ielem) const;

    int getNumEntities( const int inputDim ) const {
      if (inputDim == dim ) {
        return nels;
      } else if ( inputDim == (dim-1) ) {
        return con_FacetCell.nels_oDim;
      } else if ( inputDim == 0 ) {
        return nnodes;
      } else {
        throw std::invalid_argument("Not ready yet.");
      }
    }
    void buildAABBTree();
    void updateAABBTree();
    vector<int> findOwnerElements( const Eigen::Vector3d &point ) const;
    vector<int> findCollidingElements( const myOBB &obb ) const;
    vector<int> findCollidingElements( const Eigen::Vector3d &center, const double R) const;
};
void inPlaneRotate( Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> &points, 
    const Eigen::Vector3d &center, double angle );

}
#endif
