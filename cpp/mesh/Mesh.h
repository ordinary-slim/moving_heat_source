#ifndef MESH
#define MESH
#include <iostream>
#include <vector>
#include <Eigen/Core>
#include "Connectivity.h"
#include "Element.h"
#include "RefElement.h"
#include "ElementTypes.h"
#include "AABB.h"
#include "../../external/pybind11/include/pybind11/pybind11.h"
#include "MeshTag.h"

using namespace std;

namespace mesh
{
class Mesh {
  public:
    // Copy constructor
    Mesh(const Mesh& mesh) = default;
    // Construct from python dict with points, cells and cell_type
    Mesh(const pybind11::dict &input);
    int dim;
    int nels, nnodes, nnodes_per_el;
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>
      pos, posFRF;// node positions in xi and x
    Eigen::Vector3d shiftFRF = Eigen::Vector3d::Zero();//pos + shift = posFRF
    Eigen::Vector3d speedFRF = Eigen::Vector3d::Zero();//domain speed with respect to FRF
    Connectivity  con_CellPoint;
    Connectivity  con_PointCell;
    Connectivity  con_CellCell;
    Connectivity  con_FacetPoint;
    Connectivity  con_FacetCell;
    Connectivity  con_CellFacet;
    vector<ElementType> elementTypes;
    ReferenceElement refCellEl;
    ReferenceElement refFacetEl;
    vector<AABB> elementAABBs;

    Element getEntity(int ient, const Connectivity &connectivity, const ReferenceElement &refEl ) const;
    Element getElement(int ielem) const;

    void setSpeedFRF(Eigen::Vector3d inputSpeedFRF){
      speedFRF = inputSpeedFRF;
    }
    int getNumEntities( const int inputDim ) const {
      if (inputDim == dim ) {
        return nels;
      } else if ( inputDim == (dim-1) ) {
        return con_FacetCell.nels_oDim;
      } else if ( inputDim == 0 ) {
        return nnodes;
      } else {
        cout << "Not ready yet!" << endl;
        exit(-1);
      }
    }
    void setAABBs();
    vector<int> findOwnerElement( const Eigen::Vector3d &point );
};
MeshTag<int> mark( const Mesh &mesh, int dim = 0, const std::vector<int> &indices = std::vector<int>() );
}
#endif
