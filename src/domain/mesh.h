#ifndef MESH
#include <iostream>
#include <vector>
#include <Eigen/Core>
#include "connectivity.h"
#include "element.h"
#include "refElement.h"
#include "elementTypes.h"
#include "../../external/pybind11/include/pybind11/pybind11.h"

using namespace std;

class Mesh {
  public:
    int dim;
    int nels, nnodes, nnodes_per_el;
    Eigen::MatrixX3d pos, pos_noAdv;// node positions
    Connectivity  con_CellPoint;
    vector<ElementType> elementTypes;
    vector<int> activeElements;
    vector<int> activeNodes;
    bool hasInactive = false;
    refElement refEl;

    
    Element getElement(int ielem) {
      Element e;
      e.setElementType( refEl );

      // ALLOCATIONS
      e.pos.resize( e.nnodes, 3 );
      e.gpos.resize( e.nnodes, 3 );
      // BaseFun
      e.BaseGpVals.resize( e.nnodes );
      // Quadrature weights
      e.gpweight.resize( e.nnodes );
      // GradBaseFun
      e.GradBaseGpVals.resize( e.nnodes );
      for (int inode = 0; inode < e.nnodes; inode++) {
        e.GradBaseGpVals[inode].resize( e.ngpoints );
        for (int jgp = 0; jgp < e.ngpoints; jgp++) {
          e.GradBaseGpVals[inode][jgp].setZero();
        }
      }

      // set connectivity
      e.con = con_CellPoint.con.row(ielem);
      // set pos
      for (int inode=0; inode < e.nnodes; inode++) {
        e.pos.row(inode) = pos.row(e.con(inode));
      }

      // COMPUTATIONS
      e.computeNodalValues_Base();//COMMON BETWEEN ELS
      e.computeNodalValues_GradBase();//UNCOMMON
      return e;
    }

    void setActiveElements(vector<int> inputActiveElements ) {
      activeElements = inputActiveElements;
      //Update activeNodes
      fill( activeNodes.begin(), activeNodes.end(), 0 );
      for (int ielem = 0; ielem < nels; ielem++) {
        if (activeElements[ielem] == 1) {
          Eigen::VectorXi locCon = con_CellPoint.con.row( ielem );
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

    void generate1DMesh( double A, double B, int numberOfEls );

    // DEBUG
    void print() {
      cout << "Nodal pos:" << endl;
      for (int i = 0; i<nnodes; i++) {
        cout << pos.row(i) << endl;
      }
      cout << "connectivity:" << endl;
      for (int i = 0; i<nels; i++) {
        cout << con_CellPoint.con.row(i) << endl;
      }
    }
    void initializeMesh(pybind11::dict &input);
};
#define MESH
#endif
