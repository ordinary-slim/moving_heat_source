#include <iostream>
#include <vector>
#include "line.h"
#include "mesh.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
using namespace std;

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;

int main() {
  double L, R, h, m_ij;
  int nels;

  L = 1;
  nels = 4;

  Mesh mesh;
  mesh.initialize1DMesh( 0.0, L, nels);

  // Assembly
  Line l;
  SpMat M(mesh.nnodes, mesh.nnodes); // mass mat
  vector<T> M_coeffs;
  M_coeffs.reserve( 3*mesh.nnodes );
  SpMat K(mesh.nnodes, mesh.nnodes); // stiffness mat
  vector<T> K_coeffs;
  K_coeffs.reserve( 3*mesh.nnodes );
  for (int ielem = 0; ielem < mesh.nels; ielem++ ) {
    l = mesh.getElement( ielem );
    for (int inode = 0; inode < l.nnodes; inode++) {
      for (int jnode = 0; jnode < l.nnodes; jnode++) {
        for (int igp = 0; igp < l.nnodes; igp++) {
          m_ij = 0;
          m_ij += l.baseFunGpVals[inode][igp]*l.baseFunGpVals[jnode][igp]*l.vol;
          M_coeffs.push_back( T( l.con[inode], l.con[jnode], m_ij ) );
        }
      }
    }
  }

  M.setFromTriplets( M_coeffs.begin(), M_coeffs.end() );
  cout << "M= " << M << endl;

}
