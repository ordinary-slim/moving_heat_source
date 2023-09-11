#include "../mesh/Domain.h"
#include "Form.h"

typedef Eigen::Triplet<double> T;

void mesh::Domain::computeMassMatrix() {
  massMat.resize(mesh->nnodes, mesh->nnodes);
  vector<T> massCoeffs;
  massCoeffs.reserve( 3*mesh->nnodes );

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
      int inodeGlobal =  (*e.con)[inode] ;
      for (int jnode = 0; jnode < e.nnodes; ++jnode) {
        int jnodeGlobal = (*e.con)[jnode];
        massCoeffs.push_back( T(
              inodeGlobal,
              jnodeGlobal,
              mass_loc(inode, jnode) ) );
      }
    }
  }

  // TODO: Remove this
  // Fill mass matrix: Add 1s at inactive nodes
  vector<int> inactiveNodes = activeNodes.filterIndices( [](int activeness){ return not(activeness); });
  for (int inode : inactiveNodes) {
    massCoeffs.push_back( Eigen::Triplet<double>(inode, inode, 1) );
  }

  massMat.setFromTriplets( massCoeffs.begin(), massCoeffs.end() );
}
