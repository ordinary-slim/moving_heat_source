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
  double ip;

  std::vector<BilinearForm*> bilinearForms;
  std::vector<LinearForm*> linearForms;

  MassForm massForm = MassForm( this );
  DiffusionForm diffusionForm = DiffusionForm( this );
  AdvectionForm advectionForm = AdvectionForm( this );
  bilinearForms.push_back( &diffusionForm );
  if (isAdvection) {
    bilinearForms.push_back( &advectionForm );
  }
  /*
  if (isStabilized) {
    if (isAdvection && advectionSpeed.norm() > 1e-9) {

      ASSSBilinearForm asssLhs = ASSSBilinearForm( this );
      ASSSLinearForm asssRhs = ASSSLinearForm( this );

      bilinearForms.push_back( &asssLhs );
      linearForms.push_back( &asssRhs );
    }
  }
  */
  // Source term
  mhs.pulse.setZero();
  SourceForm sourceForm = SourceForm( this );
  linearForms.push_back( &sourceForm );

  // matrices assembly
  M.setZero();

  mesh::Element e;

  vector<int> activeElementsIndices = domain.activeElements.getTrueIndices();
  for (int ielem : activeElementsIndices) {

    e = domain.getElement( ielem );

    Eigen::MatrixXd mass_loc = Eigen::MatrixXd::Zero( e.nnodes, e.nnodes );
    Eigen::MatrixXd lhs_loc = Eigen::MatrixXd::Zero( e.nnodes, e.nnodes );
    Eigen::VectorXd rhs_loc = Eigen::VectorXd::Zero( e.nnodes );

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
        for (int jnode = 0; jnode < e.nnodes; jnode++) {
          mass_loc(inode, jnode) += massForm.contribute( igp, inode, jnode, &e );
          for (auto biForm : bilinearForms) {
            lhs_loc(inode, jnode) += biForm->contribute( igp, inode, jnode, &e );
          }
        }
      }
    }

    for (int inode = 0; inode < e.nnodes; ++inode) {
      mhs.pulse[(*e.con)[inode]] += rhs_loc(inode);
      for (int jnode = 0; jnode < e.nnodes; ++jnode) {
        massCoeffs.push_back( T( (*e.con)[inode], (*e.con)[jnode], mass_loc(inode, jnode) ) );
        lhsCoeffs.push_back( T( (*e.con)[inode], (*e.con)[jnode], lhs_loc(inode, jnode) ) );
      }
    }
  }

  M.setFromTriplets( massCoeffs.begin(), massCoeffs.end() );
  rhs += mhs.pulse;
}
