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
    vector<int> con;
    double vol;
    vector<vector<double>> baseFunGpVals;

    void setClosedIntegration(){
      gpos.reserve( nnodes );
      baseFunGpVals.resize( nnodes );
      // gauss point positions
      for (int inode = 0; inode < nnodes; inode++) {
        gpos[inode] = pos[inode];
        baseFunGpVals[inode].resize( nnodes );
        fill( baseFunGpVals[inode].begin(), baseFunGpVals[inode].end(), 0.0);
        baseFunGpVals[inode][inode] = 1.0;
      }
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
