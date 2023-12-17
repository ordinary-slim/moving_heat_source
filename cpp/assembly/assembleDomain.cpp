#include "../linearAlgebra/LinearSystem.h"
#include "../Problem.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <memory>
#include "Form.h"

typedef Eigen::Triplet<double> T;

void Problem::assembleDomain() {

  std::vector<std::unique_ptr<BilinearForm>> bilinearForms;
  std::vector<std::unique_ptr<BilinearForm>> timeDerivForms;
  std::vector<std::unique_ptr<LinearForm>> linearForms;

  bilinearForms.push_back(std::make_unique<DiffusionForm>(this));
  if (advectionSpeed.norm() > 1e-9) {
    bilinearForms.push_back(std::make_unique<AdvectionForm>(this));
    switch (stabilizationScheme) {
      case 1:
        timeDerivForms.push_back(std::make_unique<SUPGTimeBilinearForm>(this));
        bilinearForms.push_back(std::make_unique<SUPGBilinearForm>(this));
        linearForms.push_back(std::make_unique<SUPGLinearForm>(this));
        break;
    }
  }
  timeDerivForms.push_back(std::make_unique<TimeMassForm>(this));

  // Clean-up previous DSs
  timeDerivMat.resize( domain.mesh->nnodes, domain.mesh->nnodes);
  timeDerivCoeffs.clear();
  timeDerivCoeffs.reserve( 3*domain.mesh->nnodes );

  if (mhs->type == heat::lumped) {
    linearForms.push_back( std::make_unique<LumpedSourceForm>( this ) );
  } else {
    linearForms.push_back( std::make_unique<SourceForm>( this ) );
  }
  LinearForm* sourceForm = linearForms[linearForms.size()-1].get();

  mesh::Element e;

  vector<int> activeElementsIndices = domain.activeElements.getIndices();
  for (int ielem : activeElementsIndices) {

    e = domain.getElement( ielem );

    Eigen::MatrixXd lhs_loc = Eigen::MatrixXd::Zero( e.nnodes, e.nnodes );
    Eigen::MatrixXd td_loc = Eigen::MatrixXd::Zero( e.nnodes, e.nnodes );
    Eigen::VectorXd rhs_loc = Eigen::VectorXd::Zero( e.nnodes );
    Eigen::VectorXd pulse_loc = Eigen::VectorXd::Zero( e.nnodes );

    for (auto& lform : linearForms) {
      lform->preGauss( &e );
    }
    for (auto& bform : bilinearForms) {
      bform->preGauss( &e );
    }
    for (auto& tdbform : timeDerivForms ) {
      tdbform->preGauss( &e );
    }

    for (int igp = 0; igp < e.ngpoints; igp++) {
      for (auto& lform : linearForms) {
        lform->inGauss( igp, &e );
      }
      for (auto& bform : bilinearForms) {
        bform->inGauss( igp, &e );
      }
      for (auto& tdbform : timeDerivForms ) {
        tdbform->inGauss( igp, &e );
      }

      for (int inode = 0; inode < e.nnodes; inode++) {
        for (auto& lform : linearForms) {
          rhs_loc(inode) += lform->contribute( igp, inode, &e );
        }
        pulse_loc(inode) += sourceForm->contribute( igp, inode, &e );
        for (int jnode = 0; jnode < e.nnodes; jnode++) {
          for (auto& biForm : bilinearForms) {
            lhs_loc(inode, jnode) += biForm->contribute( igp, inode, jnode, &e );
          }
          for (auto& tdbform : timeDerivForms ) {
            td_loc(inode, jnode) += tdbform->contribute( igp, inode, jnode, &e );
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
        int jnodeGlobal = (*e.con)[jnode] ;;
        int jnodeDof = dofNumbering[ jnodeGlobal ];
        if ( jnodeDof < 0 ) {
          // Assemble into RHS
          ls->rhs[inodeDof] += -lhs_loc( inode, jnode ) * unknown.values[ jnodeGlobal ];
        } else {
          // Assemble into LHS
          ls->lhsCoeffs.push_back( T(
                inodeDof,
                jnodeDof,
                lhs_loc(inode, jnode) ) );
        }
        timeDerivCoeffs.push_back( T(
              inodeGlobal,
              jnodeGlobal,
              td_loc( inode, jnode ) ) );
      }
    }

    // TODO: Move / Remove this somewhere else
    // Assemble pulse
    for (int inode = 0; inode < e.nnodes; ++inode) {
      int inodeGlobal =  (*e.con)[inode] ;
      mhs->pulse[inodeGlobal] += pulse_loc(inode);
    }
  }
  
  timeDerivMat.setFromTriplets( timeDerivCoeffs.begin(),
      timeDerivCoeffs.end() );
}
