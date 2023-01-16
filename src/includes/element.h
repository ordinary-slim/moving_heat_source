#ifndef ELEMENT
#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

class Element {
  public:
    int nnodes;
    vector<double> pos, gpos;
    vector<double> gpweight;
    vector<int> con;
    double vol;
    int dimension, elementType;
    vector<vector<double>> baseFunGpVals;
    vector<vector<vector<double>>> baseFunGradGpVals;

    void setClosedIntegration(){
      //COMMON
      // gauss point positions
      for (int inode = 0; inode < nnodes; inode++) {
        gpos[inode] = pos[inode];
        baseFunGpVals[inode].resize( nnodes );
        fill( baseFunGpVals[inode].begin(), baseFunGpVals[inode].end(), 0.0);
        baseFunGpVals[inode][inode] = 1.0;
      }
      fill( gpweight.begin(), gpweight.end(), 1.0 / nnodes );

      //UNCOMMON
      setBaseFunGradGpVals();
    }

    void allocateBaseFunGradGpVals() {
      // Refactoring needed here!
      // allocate. this is common between elements
      baseFunGradGpVals.resize( nnodes );
      for (int igp = 0; igp < nnodes; igp++) {
        baseFunGradGpVals[igp].resize( nnodes );
        for (int jgp = 0; jgp < nnodes; jgp++) {
          baseFunGradGpVals[igp][jgp].resize( dimension );
        }
      }
    }

    void setBaseFunGradGpVals() {
      switch (elementType) {
        case 0: {//P0-line
          for (int jgp=0; jgp<nnodes; jgp++) {
            // compute vol
            vol = pos[ 1 ] - pos[ 0 ];
            baseFunGradGpVals[0][jgp][0] = -1.0 / vol ;
            baseFunGradGpVals[1][jgp][0] = +1.0 / vol ;
          }
          break;}
        default: {
          break;}
      }
    }

    // DEBUGGING FUNCS
    void printBaseFunGradGpVals() {
      for (int igp = 0; igp < nnodes; igp++) {
        for (int jgp = 0; jgp < nnodes; jgp++) {
          cout << "(";
          for (int idim = 0; idim < dimension; idim++) {
            printf("%.3f, ", baseFunGradGpVals[igp][jgp][idim]);
          }
          cout << ")\t\t";
        }
        cout << endl;
      }
      cout << "-------------\n";
    }

    void print() {
      for (int i = 0; i < nnodes; i++) {
        cout << pos[i] << ", ";
      }
      cout << endl;
    }
};
#define ELEMENT
#endif
