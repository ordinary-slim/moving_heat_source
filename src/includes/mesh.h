#ifndef MESH
#include <iostream>
#include <vector>
#include <Eigen/Core>
#include "element.h"

using namespace std;

class Mesh {
  public:
    int dim;
    int nels, nnodes, nnodes_per_el;
    Eigen::MatrixX3d pos;// node positions
    Eigen::MatrixXi  con;// connectivy
    vector<int> elementTypes;
    
    Element getElement(int ielem) {
      Element e;
      e.setElementType( elementTypes[ielem] );

      // ALLOCATIONS
      e.pos.resize( e.nnodes, 3 );
      e.gpos.resize( e.nnodes, 3 );
      // BaseFun
      e.BaseGpVals.resize( e.nnodes );
      // Quadrature weights
      e.gpweight.resize( e.nnodes );
      // GradBaseFun
      e.GradBaseGpVals.resize( e.nnodes );
      for (int igp = 0; igp < e.nnodes; igp++) {
        e.GradBaseGpVals[igp].resize( e.nnodes );
        for (int jgp = 0; jgp < e.nnodes; jgp++) {
          e.GradBaseGpVals[igp][jgp].setZero();
        }
      }
      // reference to current mapping
      e.ref2Local.resize( e.dim, e.dim );

      // set connectivity
      e.con = con.row(ielem);
      // set pos
      for (int inode=0; inode < e.nnodes; inode++) {
        e.pos.row(inode) = pos.row(e.con[inode]);
      }

      // COMPUTATIONS
      e.computeNodalValues_Base();//COMMON BETWEEN ELS
      e.computeNodalValues_GradBase();//UNCOMMON
      return e;
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
        cout << con.row(i) << endl;
      }
    }
};
#define MESH
#endif
