#ifndef MESH
#include <iostream>
#include <vector>
#include <Eigen/Core>
#include "Connectivity.h"
#include "Element.h"
#include "RefElement.h"
#include "ElementTypes.h"
#include "AABB.h"
#include "../../external/pybind11/include/pybind11/pybind11.h"

using namespace std;

namespace mesh
{
class Mesh {
  public:
    int dim;
    int nels, nnodes, nnodes_per_el;
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>
      pos, posFRF;// node positions in xi and x
    Eigen::Vector3d shiftFRF;//pos + shift = posFRF
    Eigen::Vector3d speedFRF;//domain speed with respect to FRF
    Connectivity  con_CellPoint;
    Connectivity  con_CellCell;
    Connectivity  con_FacetPoint;
    Connectivity  con_FacetCell;
    Connectivity  con_CellFacet;
    vector<int>   boundaryFacets;
    vector<ElementType> elementTypes;
    vector<int> activeElements;
    vector<int> activeNodes;
    bool hasInactive = false;
    int ngpointsCell = -1;//number of gausspoints for cells elements
    ReferenceElement refCellEl;
    ReferenceElement refFacetEl;
    vector<AABB> elementAABBs;

    void initializeMesh(pybind11::dict &input);
    void setActiveElements(vector<int> inputActiveElements );

    Element getEntity(int ient, Connectivity &connectivity, ReferenceElement &refEl ) {
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
    
    Element getElement(int ielem) {
      return getEntity( ielem, con_CellPoint, refCellEl );
    }
    Element getBoundaryFacet(int ifacet) {
      // Assumed that ifacet is a boundary facet
      Element e = getEntity( ifacet, con_FacetPoint, refFacetEl );
      Element parentEl = getElement( con_FacetCell.con( ifacet, 0 ) );
      e.computeNormal( parentEl.getCentroid() );
      return e;
    }

    void setSpeedFRF(Eigen::Vector3d inputSpeedFRF){
      speedFRF = inputSpeedFRF;
    }
    void setAABBs();
    int findOwnerElement( Eigen::Vector3d point );

    // DEBUG
    void print() {
      cout << "Nodal pos:" << endl;
      for (int i = 0; i<nnodes; i++) {
        cout << pos.row(i) << endl;
      }
      cout << "connectivity:" << endl;
      for (int i = 0; i<nels; i++) {
        cout << con_CellPoint.getLocalCon( i ) << endl;
      }
    }
};
}
#define MESH
#endif
