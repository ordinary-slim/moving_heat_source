#include "../Problem.h"

void Problem::assembleNeumann() {

  Eigen::VectorXd neumannRhs;
  neumannRhs.resize( mesh.nnodes );
  neumannRhs.setZero();

  double n_i;
  double flux;
  int ifacet;
  double k = material["k"];

  mesh::Element e;
  for ( int i = 0; i < neumannFacets.size(); ++i) {

    ifacet = neumannFacets[i];
    flux = neumannFluxes[i];

    e = mesh.getBoundaryFacet( ifacet );

    for (int inode = 0; inode < e.nnodes; ++inode) {
      n_i = 0;
      for (int igp = 0; igp < e.ngpoints; ++igp) {
        n_i += e.gpweight[igp] * e.vol * k * e.BaseGpVals[inode][igp] * flux;
      }
      neumannRhs[e.con[inode]] += n_i;
    }
  }

  rhs += neumannRhs;
}
