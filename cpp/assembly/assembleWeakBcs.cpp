#include "../Problem.h"
#include "Form.h"
#include <Eigen/Core>
#include <vector>

typedef Eigen::Triplet<double> T;

void Problem::assembleWeakBcs() {

  ConvectionBilinearForm convectionBi = ConvectionBilinearForm( this );
  ConvectionLinearForm convectionLi = ConvectionLinearForm( this );
  NeumannLinearForm neumannLi = NeumannLinearForm( this );

  vector<BilinearForm*> bilinearForms;
  vector<LinearForm*> linearForms;

  mesh::Element e;

  // assemble
  vector<int> indicesBoundaryFacets = domain.boundaryFacets.getTrueIndices();
  for ( int ifacet : indicesBoundaryFacets ) {

    bilinearForms.clear();
    linearForms.clear();
    // Decide what to do with facet
    if (convectionFacets[ifacet]) {
      bilinearForms.push_back( &convectionBi );
      linearForms.push_back( &convectionLi );
    } else if (neumannFacets[ifacet]) {
      linearForms.push_back( &neumannLi );
    } else {
      continue;
    }
    if (not(domain.activeElements[ domain.boundaryFacetsParentEls[ifacet] ])) {
      continue;
    }

    e = domain.getBoundaryFacet( ifacet );

    Eigen::MatrixXd lhs_loc = Eigen::MatrixXd::Zero( e.nnodes, e.nnodes );
    Eigen::VectorXd rhs_loc = Eigen::VectorXd::Zero( e.nnodes );


    for (auto lform : linearForms) {
      lform->preGauss( &e );
    }
    for (auto bform : bilinearForms) {
      bform->preGauss( &e );
    }

    for (int igp = 0; igp < e.ngpoints; ++igp) {

      // In-gauss
      for (auto lform : linearForms) {
        lform->inGauss( igp, &e );
      }
      for (auto bform : bilinearForms) {
        bform->inGauss( igp, &e );
      }

      for (int inode = 0; inode < e.nnodes; ++inode) {

        // ls.rhs-contrib
        for (auto lform : linearForms ) {
          rhs_loc(inode) += lform->contribute( igp, inode, &e );
        }

        for (int jnode = 0; jnode < e.nnodes; jnode++) {

          // rhs-contrib
          for (auto bform : bilinearForms ) {
            lhs_loc(inode, jnode) += bform->contribute( igp, inode, jnode, &e );
          }

        }
      }
    }

    for (int inode = 0; inode < e.nnodes; ++inode) {
      ls.rhs[ls.dofNumbering[(*e.con)[inode]]] += rhs_loc(inode);
      for (int jnode = 0; jnode < e.nnodes; ++jnode) {
        ls.lhsCoeffs.push_back( T(
              ls.dofNumbering[(*e.con)[inode]],
              ls.dofNumbering[(*e.con)[jnode]],
              lhs_loc(inode, jnode) ) );
      }
    }
  }
}
