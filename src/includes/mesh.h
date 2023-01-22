#ifndef MESH
#include <iostream>
#include <vector>
#include "element.h"

using namespace std;

class Mesh {
  public:
    int nels, nnodes, nnodes_per_el;
    vector<Eigen::Vector3d> pos;// node positions
    vector<vector<int>> con;// connectivy
    vector<int> elementTypes;
    
    Element getElement(int ielem) {
      Element e;
      e.elementType = elementTypes[ielem];
      switch (e.elementType) {
        case 0: {//P0-line
            e.dimension = 1;
            e.nnodes = 2;
            break;
          }
        default: {
          break;}
      }

      // ALLOCATIONS
      e.pos.reserve( e.nnodes );
      e.gpos.reserve( e.nnodes );
      // BaseFun
      e.baseFunGpVals.resize( nnodes );
      // Quadrature weights
      e.gpweight.resize( nnodes );
      // GradBaseFun
      e.baseFunGradGpVals.resize( e.nnodes );
      for (int igp = 0; igp < e.nnodes; igp++) {
        e.baseFunGradGpVals[igp].resize( e.nnodes );
        for (int jgp = 0; jgp < e.nnodes; jgp++) {
          e.baseFunGradGpVals[igp][jgp].resize( e.dimension );
        }
      }

      // set pos
      for (int inode=0; inode < e.nnodes; inode++) {
        e.pos[inode] = pos[ con[ielem][inode] ];
      }
      // set connectivity
      e.con = con[ielem];

      // COMPUTATIONS
      e.computeNodalValues_Base();//COMMON BETWEEN ELS
      e.computeNodalValues_GradBase();//UNCOMMON
      return e;
    }

    void initialize1DMesh( double A, double B, int numberOfEls );

    // DEBUG
    void print() {
      cout << "Nodal pos:" << endl;
      for (int i = 0; i<nnodes; i++) {
        cout << pos[i] << endl;
      }
      cout << "connectivity:" << endl;
      for (int i = 0; i<nels; i++) {
        for (int j = 0; j < 2; j++) {
          cout << con[i][j] << ", ";
        }
        cout << endl;
      }
    }
};
#define MESH
#endif
