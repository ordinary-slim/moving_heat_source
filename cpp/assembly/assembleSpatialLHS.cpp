#include "../Problem.h"
#include <Eigen/Core>
#include <Eigen/Sparse>

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;

void Problem::assembleSpatialLHS() {
  // numerical params
  double m_ij, k_ij, a_ij, ip;

  vector<T> K_coeffs;
  vector<T> A_coeffs;

  K_coeffs.reserve( 3*domain.mesh->nnodes );
  A_coeffs.reserve( 3*domain.mesh->nnodes );

  // matrices assembly
  M.setZero();

  mesh::Element e;
  double rho = material["rho"];
  double cp = material["cp"];
  double k = material["k"];

  vector<int> activeElementsIndices = domain.activeElements.getTrueIndices();
  for (int ielem : activeElementsIndices) {

    e = domain.getElement( ielem );

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
        }
        m_ij *= rho * cp ;
        k_ij *= k ;
        a_ij *= rho * cp;

        massCoeffs.push_back( T( (*e.con)[inode], (*e.con)[jnode], m_ij ) );
        K_coeffs.push_back( T( (*e.con)[inode], (*e.con)[jnode], k_ij ) );
        A_coeffs.push_back( T( (*e.con)[inode], (*e.con)[jnode], a_ij ) );
      }
    }
  }

  lhsCoeffs.insert(lhsCoeffs.end(), K_coeffs.begin(), K_coeffs.end());
  if (isAdvection) lhsCoeffs.insert(lhsCoeffs.end(), A_coeffs.begin(), A_coeffs.end());

  M.setFromTriplets( massCoeffs.begin(), massCoeffs.end() );
}
