#include "../Problem.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <algorithm>

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;

void Problem::assembleStabilization() {
  SpMat asss_LHS;
  Eigen::VectorXd asss_RHS;
  asss_LHS.resize( domain.mesh->nnodes, domain.mesh->nnodes );
  asss_RHS.resize( domain.mesh->nnodes );
  asss_LHS.setZero();
  asss_RHS.setZero();

  vector<T> ASSS_lhsCoeffs;
  ASSS_lhsCoeffs.reserve( 3*domain.mesh->nnodes );

  double norm_advectionSpeed = advectionSpeed.norm();

  if (!isAdvection || norm_advectionSpeed < 1e-9) {//nothing to stabilize
    return;
  }

  double SCA = 2, SCD = 4;//stabilization cte advection / diffusion
  // numerical params
  double asss_lhs_ij, asss_rhs_i;
  double tau, advectionEstimate, diffusionEstimate;//stabilization parameter
  double h;// size of element in advection direction
  double rho = material["rho"];
  double cp = material["cp"];
  double k = material["k"];
  double f_xgp;
  Eigen::Vector3d x_gp;

  mesh::Element e;
  vector<int> activeElementsIndices = domain.activeElements.getTrueIndices();

  for (int ielem : activeElementsIndices ) {

    e = domain.getElement( ielem );

    e = domain.getElement( ielem );

    //Compute tau
    h = e.getSizeAlongVector( advectionSpeed );
    advectionEstimate = h / SCA / (rho*cp*norm_advectionSpeed);
    tau = advectionEstimate;
    if (k != 0) {
      diffusionEstimate = pow(h, 2) / (SCD * k);
      tau = 1 / ( 1/advectionEstimate + 1/diffusionEstimate );
    } else {
      tau = advectionEstimate;
    }

    // Loop LHS
    for (int inode = 0; inode < e.nnodes; inode++) {
      for (int jnode = 0; jnode < e.nnodes; jnode++) {
        asss_lhs_ij = 0;
        for (int igp = 0; igp < e.ngpoints; igp++) {
          asss_lhs_ij += (e.gpweight[igp] * e.vol)* rho * cp * tau *
            e.GradBaseGpVals[inode][igp].dot( advectionSpeed ) *
            e.GradBaseGpVals[jnode][igp].dot( advectionSpeed );
        }
        ASSS_lhsCoeffs.push_back( T( (*e.con)[inode], (*e.con)[jnode], asss_lhs_ij ) );
      }
    }

    // Loop RHS
    for (int inode = 0; inode < e.nnodes; inode++) {
      asss_rhs_i = 0;
      for (int igp = 0; igp < e.ngpoints; igp++) {
        x_gp = e.gpos.row( igp );
        f_xgp = mhs.powerDensity(x_gp, time, mhs.currentPosition, mhs.power, mhs.efficiency, mhs.radius);
        asss_rhs_i += (e.gpweight[igp] * e.vol) * tau * f_xgp *
          e.GradBaseGpVals[inode][igp].dot( advectionSpeed );;
      }
      asss_RHS[(*e.con)[inode]] += asss_rhs_i;
    }

  }

  //asss_LHS.setFromTriplets( LHS_coeffs.begin(), LHS_coeffs.end() );

  lhsCoeffs.insert( lhsCoeffs.end(), ASSS_lhsCoeffs.begin(), ASSS_lhsCoeffs.end());
  rhs += asss_RHS;
}
