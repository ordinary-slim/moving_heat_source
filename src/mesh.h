#ifndef MESH
#include <iostream>
#include <vector>
#include "line.h"

using namespace std;

class Mesh {
  public:
    int nels, nnodes, nnodes_per_el;
    vector<double> pos;// node positions
    vector<vector<int>> con;// connectivy
    
    void initialize1DMesh( double A, double B, int numberOfEls ){
      nels = numberOfEls;
      nnodes = nels + 1;
      double L = (B - A);
      double h = L / nels;

      con.resize( nels );
      pos.resize( nnodes );

      for (int i = 0; i < nnodes; i++) {
        pos[i] = i*h;
      }

      // build connectivity; inneficient but clear
      nnodes_per_el = 2;
      for (int i = 0; i < nels; i++) {
        con[i].push_back( i );
        con[i].push_back( i+1 );
      }
    }

    Line getElement(int ielem) {
      Line l;
      l.pos[0] = pos[ con[ielem][0] ];
      l.pos[1] = pos[ con[ielem][1] ];
      l.vol = pos[ con[ielem][1] ] - pos[ con[ielem][0] ];
      l.con = con[ielem];
      l.setClosedIntegration();
      return l;
    }

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
