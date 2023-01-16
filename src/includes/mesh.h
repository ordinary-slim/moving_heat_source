#ifndef MESH
#include <iostream>
#include <vector>
#include "element.h"

using namespace std;

class Mesh {
  public:
    int nels, nnodes, nnodes_per_el;
    vector<double> pos;// node positions
    vector<vector<int>> con;// connectivy
    vector<int> elementTypes;
    
    void initialize1DMesh( double A, double B, int numberOfEls );

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

      // allocs
      e.pos.reserve( e.nnodes );
      e.gpos.reserve( e.nnodes );
      e.baseFunGpVals.resize( nnodes );
      e.gpweight.resize( nnodes );
      e.allocateBaseFunGradGpVals();

      // set pos
      for (int inode=0; inode < e.nnodes; inode++) {
        e.pos[inode] = pos[ con[ielem][inode] ];
      }
      // set connectivity
      e.con = con[ielem];

      e.setClosedIntegration();//derivatives etc
      return e;
    }

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
