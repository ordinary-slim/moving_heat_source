#include "../Problem.h"

void Problem::assembleNeumann() {

  Eigen::VectorXd neumannRhs;
  neumannRhs.resize( domain.mesh->nnodes );
  neumannRhs.setZero();

  double n_i;
  double k = material["k"];
  double normalDerivative;

  mesh::Element e;

  vector<int> indicesNeumanFacets = neumannFacets.getTrueIndices();
  for ( int ifacet : indicesNeumanFacets ) {
    e = domain.getBoundaryFacet( ifacet );

    for (int inode = 0; inode < e.nnodes; ++inode) {
      n_i = 0;
      for (int igp = 0; igp < e.ngpoints; ++igp) {
        normalDerivative = - neumannFluxes[ifacet][igp] / k;
        n_i += e.gpweight[igp] * e.vol * e.BaseGpVals[inode][igp] * normalDerivative;
      }
      neumannRhs[(*e.con)[inode]] += n_i;
    }
  }

  rhs += neumannRhs;
}
