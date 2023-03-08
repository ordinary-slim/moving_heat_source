#include "Problem.h"
#include <Eigen/Core>
#include <Eigen/Sparse>

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;

void Problem::assemble() {
  // numerical params
  double m_ij, k_ij, a_ij, ip;

  if (!isAssembled) {
    // initialize data structures
    M.resize(mesh.nnodes, mesh.nnodes); // mass mat
    K.resize(mesh.nnodes, mesh.nnodes); // stiffness mat
    A.resize(mesh.nnodes, mesh.nnodes); // advection mat
    I.resize(mesh.nnodes, mesh.nnodes); // inactive nodes
    lhs.resize( mesh.nnodes, mesh.nnodes );
    rhs.resize( mesh.nnodes );
    pulse.resize( mesh.nnodes ); // source term

    vector<T> M_coeffs;
    M_coeffs.reserve( 3*mesh.nnodes );
    vector<T> K_coeffs;
    K_coeffs.reserve( 3*mesh.nnodes );
    vector<T> A_coeffs;
    A_coeffs.reserve( 3*mesh.nnodes );

    // matrices assembly
    M.setZero();
    K.setZero();
    A.setZero();
    I.setZero();

    mesh::Element e;
    double rho = material["rho"];
    double cp = material["cp"];
    double k = material["k"];
    for (int ielem = 0; ielem < mesh.nels; ielem++ ) {
      if (mesh.activeElements[ielem]==0){
        continue;
      }
      e = mesh.getElement( ielem );
      for (int inode = 0; inode < e.nnodes; inode++) {
        for (int jnode = 0; jnode < e.nnodes; jnode++) {
          m_ij = 0;
          k_ij = 0;
          a_ij = 0;
          for (int igp = 0; igp < e.ngpoints; igp++) {
            // mass matrix
            m_ij += e.gpweight[igp] * e.BaseGpVals[inode][igp]*e.BaseGpVals[jnode][igp]*e.vol;
            // stiffness matrix
            ip = e.GradBaseGpVals[inode][igp].dot( e.GradBaseGpVals[jnode][igp] );
            k_ij += e.gpweight[igp] * ip * e.vol;
            // advection matrix
            ip = e.GradBaseGpVals[jnode][igp].dot(advectionSpeed);
            a_ij += e.gpweight[igp] * (ip * e.BaseGpVals[inode][igp]) * e.vol;
            // stabilization matrix
            // TODO: implement computation of s_ij @ gauss point
          }
          m_ij *= rho * cp ;
          k_ij *= k ;
          a_ij *= rho * cp;
          //TODO: s_ij *= blablabla;

          // TODO: Implement Dirichlet BCs
          // lookup i or j belong to fixed nodes
          //
          // if belong
          // send it to rhs multiplying DIRICHLET


          M_coeffs.push_back( T( e.con[inode], e.con[jnode], m_ij ) );
          K_coeffs.push_back( T( e.con[inode], e.con[jnode], k_ij ) );
          A_coeffs.push_back( T( e.con[inode], e.con[jnode], a_ij ) );
        }
      }
    }
    M.setFromTriplets( M_coeffs.begin(), M_coeffs.end() );
    K.setFromTriplets( K_coeffs.begin(), K_coeffs.end() );
    A.setFromTriplets( A_coeffs.begin(), A_coeffs.end() );

    isAssembled = true;
  }
}
