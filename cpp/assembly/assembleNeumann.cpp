#include "../Problem.h"

void Problem::assembleNeumann() {

  Eigen::VectorXd neumannRhs;
  neumannRhs.resize( domain.mesh->nnodes );
  neumannRhs.setZero();

  double n_i;
  double normalDerivative;
  int ifacet;
  double k = material["k"];

  mesh::Element e;
  for ( int i = 0; i < neumannFacets.size(); ++i) {

    ifacet = neumannFacets[i];

    e = domain.getBoundaryFacet( ifacet );

    for (int inode = 0; inode < e.nnodes; ++inode) {
      n_i = 0;
      for (int igp = 0; igp < e.ngpoints; ++igp) {
        normalDerivative = - neumannFluxes[i][igp] / k;
        n_i += e.gpweight[igp] * e.vol * e.BaseGpVals[inode][igp] * normalDerivative;
      }
      neumannRhs[(*e.con)[inode]] += n_i;
    }
  }

  rhs += neumannRhs;
}
