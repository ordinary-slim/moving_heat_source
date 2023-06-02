#include "../Problem.h"

typedef Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> Dense3ColMat;
typedef Eigen::Triplet<double> T;

void Problem::assembleNeumannGamma(const Problem &pExt) {
  vector<int> gammaIndices = gammaFacets.getIndices();

  for (int ifacet : gammaIndices) {
    mesh::Element facet = domain.getBoundaryFacet( ifacet );

    for (int igp = 0; igp < facet.ngpoints; ++igp) {

      Eigen::Vector3d xgp = facet.gpos.row( igp );
      int idx_el_ext = pExt.domain.findOwnerElement( xgp );
      mesh::Element e_ext = pExt.domain.getElement( idx_el_ext );

      Dense3ColMat gradShaFuns = e_ext.evaluateGradShaFuns( xgp );

      Eigen::MatrixXd lhs_loc = Eigen::MatrixXd::Zero( facet.nnodes, e_ext.nnodes );

      for (int jnode = 0; jnode < e_ext.nnodes; jnode++) {

        double gradj_n = gradShaFuns.row( jnode ).dot( facet.normal );

        for (int inode = 0; inode < facet.nnodes; ++inode) {

          lhs_loc(inode, jnode) += -pExt.conductivity * 
                          facet.BaseGpVals[inode][igp] * gradj_n *
                          facet.gpweight[igp] * facet.vol;
        }
      }

      for (int inode = 0; inode < facet.nnodes; ++inode) {

        int inodeGlobal =  (*facet.con)[inode];
        int inodeDof = dofNumbering[ inodeGlobal ];
        if ( inodeDof < 0 ) { continue; }// if forced node, keep going

        for (int jnode = 0; jnode < e_ext.nnodes; jnode++) {

          int jnodeGlobal_ext =  (*e_ext.con)[jnode];
          int jnodeDof = pExt.dofNumbering[ jnodeGlobal_ext ];

          if ( jnodeDof < 0 ) {
            // Assemble into RHS
            ls->rhs[inodeDof] += - lhs_loc( inode, jnode ) * pExt.unknown.values[ jnodeGlobal_ext ];
          } else {
            // Assemble into RHS
            ls->lhsCoeffs.push_back( T(
                  inodeDof,
                  jnodeDof,
                  lhs_loc(inode, jnode) ) );
          }
        }
      }

    }

  }

}

void Problem::assembleDirichletGamma(const Problem &pExt) {
}
