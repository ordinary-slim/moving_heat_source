#ifndef ELEMENT
#include <iostream>
#include <vector>
#include <Eigen/Core>
#include <algorithm>
#include "refElement.h"
using namespace std;

class Element {
  public:
    int nnodes, ngpoints;
    Eigen::MatrixX3d pos, gpos;
    Eigen::VectorXi  con;
    vector<double> gpweight;
    double vol;
    int dim, elementType;
    vector<vector<double>> BaseGpVals;
    vector<vector<Eigen::Vector3d>> GradBaseGpVals;

    refElement      *refEl;
    Eigen::Matrix3d ref2loc;// x = c + ref2loc · xi
    Eigen::Matrix3d loc2ref_T;// tranpose of loc2ref where loc2ref st xi= c'+ loc2ref · x

    void computeNodalValues_Base(){
      //COMMON
      //Closed integration
      for (int inode = 0; inode < nnodes; inode++) {
        gpos.row(inode) = pos.row(inode);
        BaseGpVals[inode].resize( nnodes );
        fill( BaseGpVals[inode].begin(), BaseGpVals[inode].end(), 0.0);
        BaseGpVals[inode][inode] = 1.0;
      }
      fill( gpweight.begin(), gpweight.end(), 1.0 / nnodes );
    }

    void computeNodalValues_GradBase() {
      /*
      switch (elementType) {
        case 0: {//P0-line
          vol = pos(1, 0) - pos(0, 0);
          for (int jgp=0; jgp<nnodes; jgp++) {
            // compute vol
            GradBaseGpVals[0][jgp][0] = -1.0 / vol ;
            GradBaseGpVals[1][jgp][0] = +1.0 / vol ;
          }
          break;}
        default: {
          break;}
      }
      */
      computeDerivatives();
    }

    void setElementType( refElement target_refEl ) {
      refEl = &target_refEl;
      elementType = refEl->elementType;
      nnodes = refEl->nnodes;
      ngpoints = refEl->ngpoints;
      dim = refEl->dim;
    }

    void computeDerivatives() {
      // Build helper matrices
      Eigen::Matrix3d X;
      // Compute
      X.setZero();
      for (int inode = 0; inode<dim; inode++) {
        X.row( inode )  = pos.row( inode+1 ) - pos.row( 0 );
      }
      X.transposeInPlace();
      // Fill missing dims
      for (int idim = dim; idim < 3; idim++) {
        X(idim, idim)  = 1.0;
      }
      ref2loc = X * refEl->XI_inverse;
      loc2ref_T = ref2loc.inverse();
      loc2ref_T.transposeInPlace();

      // Compute volume
      vol = refEl->vol * ref2loc.determinant();

      for (int inode = 0; inode < nnodes; inode++) {
        for (int igp = 0; igp < ngpoints; igp++) {
          GradBaseGpVals[inode][igp] = (loc2ref_T)*(refEl->GradBaseGpVals[inode][igp]);
        }
      }
    }


    // DEBUGGING FUNCS
    void printBaseFunGradGpVals() {
      for (int igp = 0; igp < nnodes; igp++) {
        for (int jgp = 0; jgp < nnodes; jgp++) {
          cout << "(";
          for (int idim = 0; idim < dim; idim++) {
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
        cout << pos.row(i) << ", ";
      }
      cout << endl;
    }
};
#define ELEMENT
#endif
