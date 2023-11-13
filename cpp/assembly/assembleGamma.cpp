#include "../Problem.h"

typedef Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> Dense3ColMat;
typedef Eigen::Triplet<double> T;

void Problem::assembleNeumannGamma(const Problem *pExt) {
  vector<int> gammaNeumannIndices = weakBcFacets.filterIndices( [](int tag){ return tag==3; } );

  for (int ifacet : gammaNeumannIndices) {
    mesh::Element facet = domain.getBoundaryFacet( ifacet );

    for (int igp = 0; igp < facet.ngpoints; ++igp) {

      // Position of gauss point in reference frame of external problem
      Eigen::Vector3d xgp_ext = facet.gpos.row( igp ).transpose() + domain.translationLab - pExt->domain.translationLab;
      // Find out which element of external mesh owns Gauss point
      int idx_el_ext = pExt->domain.findOwnerElements( xgp_ext );
      if (idx_el_ext == -1) { mesh::throwPointOutOfBounds( xgp_ext ); };
      mesh::Element e_ext = pExt->domain.getElement( idx_el_ext );
      const ThermalMaterial* ext_mat = &pExt->materials[ e_ext.imat ];

      Dense3ColMat gradShaFuns = e_ext.evaluateGradShaFuns( xgp_ext );

      Eigen::MatrixXd lhs_loc = Eigen::MatrixXd::Zero( facet.nnodes, e_ext.nnodes );

      for (int jnode = 0; jnode < e_ext.nnodes; jnode++) {

        double gradj_n = gradShaFuns.row( jnode ).dot( facet.normal );

        for (int inode = 0; inode < facet.nnodes; ++inode) {

          lhs_loc(inode, jnode) += -ext_mat->conductivity * 
                          facet.BaseGpVals[inode][igp] * gradj_n *
                          facet.gpweight[igp] * facet.vol;
        }
      }

      for (int inode = 0; inode < facet.nnodes; ++inode) {

        int inodeGlobal =  (*facet.con)[inode];
        int inodeDof = freeDofsNumbering[ inodeGlobal ];
        if ( inodeDof < 0 ) { continue; }// if forced node, keep going

        for (int jnode = 0; jnode < e_ext.nnodes; jnode++) {

          int jnodeGlobal_ext =  (*e_ext.con)[jnode];
          int jnodeDof = pExt->dofNumbering[ jnodeGlobal_ext ];

          if ( jnodeDof < 0 ) {
            // Assemble into RHS
            ls->rhs[inodeDof] += - lhs_loc( inode, jnode ) * pExt->unknown.values[ jnodeGlobal_ext ];
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

void Problem::assembleDirichletGamma(const Problem *pExt) {
  vector<int> gammaDirichletIndices = dirichletNodes.filterIndices( [](int tag){ return tag==2; } );
  for (int inode: gammaDirichletIndices) {
    Eigen::Vector3d xnode_ext = domain.mesh->pos.row( inode ).transpose() + domain.translationLab - pExt->domain.translationLab;

    int idx_el_ext = pExt->domain.findOwnerElements( xnode_ext ) ;
    if (idx_el_ext == -1) { mesh::throwPointOutOfBounds( xnode_ext ); };
    mesh::Element e_ext = pExt->domain.getElement( idx_el_ext );

    Eigen::VectorXd shaFuns = e_ext.evaluateShaFuns( xnode_ext );

    // Assemble
    int inodeDof = dofNumbering[inode];

    ls->lhsCoeffs.push_back( T(
          inodeDof,
          inodeDof,
          -1.0 ) ); 
    for (int jnode = 0; jnode < e_ext.nnodes; ++jnode) {
      int jnodeGlobal_ext =  (*e_ext.con)[jnode];
      int jnodeDof = pExt->dofNumbering[ jnodeGlobal_ext ];
      if ( jnodeDof < 0 ) {
        // Assemble to RHS
        ls->rhs[inodeDof] += - shaFuns[jnode] * pExt->unknown.values[ jnodeGlobal_ext ];
      } else {
        // Assemble to LHS
        ls->lhsCoeffs.push_back( T(
              inodeDof,
              jnodeDof,
              shaFuns[jnode] ) ); 
      }
    }

  }
}
