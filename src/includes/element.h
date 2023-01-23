#ifndef ELEMENT
#include <iostream>
#include <vector>
#include <Eigen/Core>
#include <algorithm>
using namespace std;

class Element {
  public:
    int nnodes;
    vector<Eigen::Vector3d> pos, gpos;
    vector<double> gpweight;
    vector<int> con;
    double vol;
    int dimension, elementType;
    vector<vector<double>> BaseGpVals;
    vector<vector<Eigen::Vector3d>> GradBaseGpVals;

    void computeNodalValues_Base(){
      //COMMON
      //Closed integration
      for (int inode = 0; inode < nnodes; inode++) {
        gpos[inode] = pos[inode];
        BaseGpVals[inode].resize( nnodes );
        fill( BaseGpVals[inode].begin(), BaseGpVals[inode].end(), 0.0);
        BaseGpVals[inode][inode] = 1.0;
      }
      fill( gpweight.begin(), gpweight.end(), 1.0 / nnodes );
    }

    void computeNodalValues_GradBase() {
      switch (elementType) {
        case 0: {//P0-line
          vol = pos[ 1 ][0] - pos[ 0 ][0];
          for (int jgp=0; jgp<nnodes; jgp++) {
            // compute vol
            GradBaseGpVals[0][jgp][0] = -1.0 / vol ;
            GradBaseGpVals[1][jgp][0] = +1.0 / vol ;
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
            printf("%.3f, ", GradBaseGpVals[igp][jgp][idim]);
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
