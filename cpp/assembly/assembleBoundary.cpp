#include "../Problem.h"
#include "Form.h"
#include <Eigen/Core>
#include <vector>
#include <memory>

typedef Eigen::Triplet<double> T;

void Problem::assembleBoundary() {

  std::vector<std::unique_ptr<BilinearForm>> bilinearForms;
  std::vector<std::unique_ptr<LinearForm>> linearForms;

  // assemble
  vector<int> indicesBcFacets = weakBcFacets.getIndices();
  for ( int ifacet : indicesBcFacets ) {

    // Support for different BCs depending on facet
    bilinearForms.clear();//resetting facet
    linearForms.clear();
    // Decide what to do with facet
    switch (weakBcFacets[ifacet]) {
      case 1:
        linearForms.push_back(std::make_unique<NeumannLinearForm>(this));
        break;
      case 2:
        bilinearForms.push_back(std::make_unique<ConvectionBilinearForm>(this));
        linearForms.push_back(std::make_unique<ConvectionLinearForm>(this));
        break;
      default:
        continue;
    }

    if (not(domain.activeElements[ domain.boundaryFacetsParentEls[ifacet] ])) {
      continue;
    }

    mesh::Element e = domain.getBoundaryFacet( ifacet );

    Eigen::MatrixXd lhs_loc = Eigen::MatrixXd::Zero( e.nnodes, e.nnodes );
    Eigen::VectorXd rhs_loc = Eigen::VectorXd::Zero( e.nnodes );

    for (auto& lform : linearForms) {
      lform->preGauss( &e );
    }
    for (auto& bform : bilinearForms) {
      bform->preGauss( &e );
    }

    for (int igp = 0; igp < e.ngpoints; ++igp) {

      // In-gauss
      for (auto& lform : linearForms) {
        lform->inGauss( igp, &e );
      }
      for (auto& bform : bilinearForms) {
        bform->inGauss( igp, &e );
      }

      for (int inode = 0; inode < e.nnodes; ++inode) {

        // rhs-contrib
        for (auto& lform : linearForms ) {
          rhs_loc(inode) += lform->contribute( igp, inode, &e );
        }

        for (int jnode = 0; jnode < e.nnodes; jnode++) {

          // rhs-contrib
          for (auto& bform : bilinearForms ) {
            lhs_loc(inode, jnode) += bform->contribute( igp, inode, jnode, &e );
          }

        }
      }
    }

    // Assemble into linear system 
    for (int inode = 0; inode < e.nnodes; ++inode) {
      int inodeGlobal =  (*e.con)[inode] ;
      int inodeDof = freeDofsNumbering[ inodeGlobal ];
      if ( inodeDof < 0 ) { continue; }// if forced node, keep going
                                       //
      // Assemble into RHS
      ls->rhs[inodeDof] += rhs_loc(inode);

      for (int jnode = 0; jnode < e.nnodes; ++jnode) {
        int jnodeGlobal = (*e.con)[jnode] ;
        int jnodeDof = dofNumbering[ jnodeGlobal ];
        if ( jnodeDof < 0 ) {
          // Assemble into RHS
          ls->rhs[inodeDof] += - lhs_loc( inode, jnode ) * unknown.values[ jnodeGlobal ];
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
