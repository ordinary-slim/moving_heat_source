#include "../Problem.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include "Form.h"

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;

double massContrib( int igp, int inode, int jnode, const mesh::Element *e, const Problem *p ) {
  return p->density * p->specificHeat *
    e->BaseGpVals[inode][igp]*e->BaseGpVals[jnode][igp]*
    e->gpweight[igp]  * e->vol;
}

void Problem::assembleSpatialPDE() {
  // numerical params

  std::vector<BilinearForm*> bilinearForms;
  std::vector<LinearForm*> linearForms;

  MassForm massForm = MassForm( this );
  DiffusionForm diffusionForm = DiffusionForm( this );
  AdvectionForm advectionForm = AdvectionForm( this );
  ASSSBilinearForm asssLhs = ASSSBilinearForm( this );
  ASSSLinearForm asssRhs = ASSSLinearForm( this );

  bilinearForms.push_back( &diffusionForm );
  if (isAdvection) {
    bilinearForms.push_back( &advectionForm );
  }
  if (isStabilized) {
    if (isAdvection && advectionSpeed.norm() > 1e-9) {

      bilinearForms.push_back( &asssLhs );
      linearForms.push_back( &asssRhs );
    }
  }

  // Source term
  mhs.pulse.setZero();
  SourceForm sourceForm = SourceForm( this );
  linearForms.push_back( &sourceForm );

  domain.massMat.setZero();

  mesh::Element e;

  vector<int> activeElementsIndices = domain.activeElements.getIndices();
  for (int ielem : activeElementsIndices) {

    e = domain.getElement( ielem );

    Eigen::MatrixXd mass_loc = Eigen::MatrixXd::Zero( e.nnodes, e.nnodes );
    Eigen::MatrixXd lhs_loc = Eigen::MatrixXd::Zero( e.nnodes, e.nnodes );
    Eigen::VectorXd rhs_loc = Eigen::VectorXd::Zero( e.nnodes );
    Eigen::VectorXd pulse_loc = Eigen::VectorXd::Zero( e.nnodes );

    for (auto lform : linearForms) {
      lform->preGauss( &e );
    }
    for (auto bform : bilinearForms) {
      bform->preGauss( &e );
    }

    for (int igp = 0; igp < e.ngpoints; igp++) {
      for (auto lform : linearForms) {
        lform->inGauss( igp, &e );
      }
      for (auto bform : bilinearForms) {
        bform->inGauss( igp, &e );
      }
      for (int inode = 0; inode < e.nnodes; inode++) {
        for (auto lform : linearForms) {
          rhs_loc(inode) += lform->contribute( igp, inode, &e );
        }
        pulse_loc(inode) += sourceForm.contribute( igp, inode, &e );
        for (int jnode = 0; jnode < e.nnodes; jnode++) {
          mass_loc(inode, jnode) += massForm.contribute( igp, inode, jnode, &e );
          for (auto biForm : bilinearForms) {
            lhs_loc(inode, jnode) += biForm->contribute( igp, inode, jnode, &e );
          }
        }
      }
    }

    // Assemble into linear system 
    for (int inode = 0; inode < e.nnodes; ++inode) {
      int inodeGlobal =  (*e.con)[inode] ;
      int inodeDof = dofNumbering[ inodeGlobal ];
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
          // Assemble into RHS
          ls->lhsCoeffs.push_back( T(
                inodeDof,
                jnodeDof,
                lhs_loc(inode, jnode) ) );
        }
      }
    }

    // TODO: Move this somewhere else
    // Assemble mass and pulse
    for (int inode = 0; inode < e.nnodes; ++inode) {
      int inodeGlobal =  (*e.con)[inode] ;
      mhs.pulse[inodeGlobal] += pulse_loc(inode);
      for (int jnode = 0; jnode < e.nnodes; ++jnode) {
        int jnodeGlobal = (*e.con)[jnode];
        domain.massCoeffs.push_back( T(
              inodeGlobal,
              jnodeGlobal,
              mass_loc(inode, jnode) ) );
      }
    }
  }

  domain.massMat.setFromTriplets( domain.massCoeffs.begin(), domain.massCoeffs.end() );
}
