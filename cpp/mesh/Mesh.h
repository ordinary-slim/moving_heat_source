#ifndef MESH
#include <iostream>
#include <vector>
#include <Eigen/Core>
#include "connectivity.h"
#include "Element.h"
#include "refElement.h"
#include "elementTypes.h"
#include "AABB.h"
#include "../../external/pybind11/include/pybind11/pybind11.h"

using namespace std;

namespace mesh
{
class Mesh {
  public:
    int dim;
    int nels, nnodes, nnodes_per_el;
    Eigen::MatrixX3d pos, posFRF;// node positions
    Connectivity  con_CellPoint;
    Connectivity  con_CellCell;
    Connectivity  con_FacetPoint;
    Connectivity  con_FacetCell;
    Connectivity  con_CellFacet;
    vector<ElementType> elementTypes;
    vector<int> activeElements;
    vector<int> activeNodes;
    bool hasInactive = false;
    ReferenceElement refCellEl;
    ReferenceElement refFacetEl;
    vector<AABB> elementAABBs;

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
    Element getFacetElement(int ifacet) {
      return getEntity( ifacet, con_FacetPoint, refFacetEl );
    }

    void setActiveElements(vector<int> inputActiveElements ) {
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

    void setAABBs();
    int findOwnerElement( Eigen::Vector3d point );
    void generate1DMesh( double A, double B, int numberOfEls );

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
    void initializeMesh(pybind11::dict &input);
};
}
#define MESH
#endif
