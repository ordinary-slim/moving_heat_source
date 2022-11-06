#ifndef ELEMENT
#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

class Element {
  public:
    int nnodes;
    vector<double> rpos, rgpos;
    vector<double> pos, gpos;
    vector<double> gpweight;
    vector<int> con;
    double vol;
    int dimension;
    vector<vector<double>> baseFunGpVals;
    vector<vector<vector<double>>> baseFunGradGpVals;

    void setClosedIntegration(){
      // common for all elements
      gpos.reserve( nnodes );
      baseFunGpVals.resize( nnodes );
      gpweight.resize( nnodes );
      // gauss point positions
      for (int inode = 0; inode < nnodes; inode++) {
        gpos[inode] = pos[inode];
        baseFunGpVals[inode].resize( nnodes );
        fill( baseFunGpVals[inode].begin(), baseFunGpVals[inode].end(), 0.0);
        baseFunGpVals[inode][inode] = 1.0;
      }

      fill( gpweight.begin(), gpweight.end(), 1.0 / nnodes );

      allocateBaseFunGradGpVals();
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

    virtual void setBaseFunGradGpVals() = 0;

    void print() {
      for (int i = 0; i < nnodes; i++) {
        cout << pos[i] << ", ";
      }
      cout << endl;
    }
};
#define ELEMENT
#endif
