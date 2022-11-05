#ifndef ELEMENT
#include <iostream>
#include <vector>
using namespace std;

class Element {
  public:
    int nnodes;
    vector<double> rpos, rgpos;
    vector<double> pos, gpos;
    vector<vector<double>> baseFunGpValues;
    double vol;

    void setGposClosed(){
      baseFunGpValues.reserve( nnodes );
      gpos.reserve( nnodes );
      // gauss point positions
      for (int inode = 0; inode < nnodes; inode++) {
        gpos[inode] = pos[inode];
        baseFunGpValues[inode].reserve( nnodes );
        for (int igp = 0; igp < nnodes; igp++) {
          baseFunGpValues[inode][igp] = 0.0;
          if ( inode==igp ) { baseFunGpValues[inode][igp] = 1.0; };
        }
      }
    }

    void print() {
      for (int i = 0; i < nnodes; i++) {
        cout << pos[i] << ", ";
      }
      cout << endl;
      // values base func at gp
      for (int inode=0; inode < nnodes; inode++) {
        for (int igp = 0; igp < nnodes; igp++) {
          printf("(%d, %d) ---- %.2f\n", inode, igp, baseFunGpValues[inode][igp]);
        }
      }
    }
};
#define ELEMENT
#endif
