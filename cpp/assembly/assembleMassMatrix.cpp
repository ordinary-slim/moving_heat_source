#include "../mesh/Domain.h"
#include "Form.h"

typedef Eigen::Triplet<double> T;

void mesh::Domain::assembleMassMatrix() {

  // Allocate Linear System
  this->ls = std::make_shared<LinearSystem>(*this);

  MassForm massForm = MassForm();

  mesh::Element e;

  vector<int> activeElementsIndices = activeElements.getIndices();
  for (int ielem : activeElementsIndices) {
    e = getElement( ielem );
    Eigen::MatrixXd mass_loc = Eigen::MatrixXd::Zero( e.nnodes, e.nnodes );
    for (int igp = 0; igp < e.ngpoints; igp++) {
      for (int inode = 0; inode < e.nnodes; inode++) {
        for (int jnode = 0; jnode < e.nnodes; jnode++) {
          mass_loc(inode, jnode) += massForm.contribute( igp, inode, jnode, &e );
        }
      }
    }

    // Assemble to global coeffs
    for (int inode = 0; inode < e.nnodes; ++inode) {
      int inodeDof =  dofNumbering[(*e.con)[inode]] ;
      if (inodeDof < 0) { continue; }// if node is inactive skip it
      for (int jnode = 0; jnode < e.nnodes; ++jnode) {
        int jnodeDof = dofNumbering[(*e.con)[jnode]];
        if (jnodeDof < 0) { continue; }// if node is inactive skip it
        ls->lhsCoeffs.push_back( T(
              inodeDof,
              jnodeDof,
              mass_loc(inode, jnode) ) );
      }
    }
  }

  ls->assemble();
  massMat = &ls->lhs;
}
