#include <iostream>
#include <vector>
#include "line.h"
#include "mesh.h"
#include <numeric>
#include <Eigen/Core>
#include <Eigen/Sparse>
using namespace std;

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;

int main() {
  // params
  double L=20.0, R;
  double x0 = 0.0, v=10, t=0.0, Tfinal=10.0;
  // numerical params
  double dt;
  double CFL;
  double m_ij, k_ij, ip;
  int nels;
  int iter=0, maxIter=100;

  nels = 4;
  Mesh mesh;
  mesh.initialize1DMesh( 0.0, L, nels);

  // set timestep
  CFL = 1.0;
  dt = CFL * (L / nels) / v;
  cout << dt << endl;


  // Assembly
  SpMat M(mesh.nnodes, mesh.nnodes); // mass mat
  vector<T> M_coeffs;
  M_coeffs.reserve( 3*mesh.nnodes );
  SpMat K(mesh.nnodes, mesh.nnodes); // stiffness mat
  vector<T> K_coeffs;
  K_coeffs.reserve( 3*mesh.nnodes );
  Eigen::VectorXd r( mesh.nnodes ); // source term

  Line l;

  while ((t < Tfinal)&&(iter < maxIter)) {
      for (int ielem = 0; ielem < mesh.nels; ielem++ ) {
        l = mesh.getElement( ielem );
        for (int inode = 0; inode < l.nnodes; inode++) {
          for (int jnode = 0; jnode < l.nnodes; jnode++) {
            for (int igp = 0; igp < l.nnodes; igp++) {
              // mass matrix
              m_ij = 0;
              m_ij += l.gpweight[igp] * l.baseFunGpVals[inode][igp]*l.baseFunGpVals[jnode][igp]*l.vol;
              M_coeffs.push_back( T( l.con[inode], l.con[jnode], m_ij ) );
              // stiffness matrix
              k_ij = 0;
              ip = inner_product(l.baseFunGradGpVals[inode][igp].begin(),
                    l.baseFunGradGpVals[inode][igp].end(),
                    l.baseFunGradGpVals[jnode][igp].begin(),
                    0.0);
              k_ij += l.gpweight[igp] * ip * l.vol;
              K_coeffs.push_back( T( l.con[inode], l.con[jnode], k_ij ) );
            }
          }
        }
      }
    ++iter;
    t += dt;
  }

  M.setFromTriplets( M_coeffs.begin(), M_coeffs.end() );
  K.setFromTriplets( K_coeffs.begin(), K_coeffs.end() );

}
