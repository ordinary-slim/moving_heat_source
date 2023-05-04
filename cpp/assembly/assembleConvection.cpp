#include "../Problem.h"

typedef Eigen::Triplet<double> T;

void Problem::assembleConvectionLHS() {

  vector<T> C_coeffs;
  C_coeffs.reserve( 1*domain.mesh->nnodes );

  double k = material["k"];
  double h = material["h"];

  mesh::Element e;
  double c_ij;

  // assemble
  vector<int> indicesConvectionFacets = convectionFacets.getTrueIndices();
  for ( int ifacet : indicesConvectionFacets ) {

    //TODO: Think about active elements!
    e = domain.getBoundaryFacet( ifacet );

    for (int inode = 0; inode < e.nnodes; ++inode) {
      for (int jnode = 0; jnode < e.nnodes; jnode++) {
        c_ij = 0.0;
        for (int igp = 0; igp < e.ngpoints; ++igp) {
          c_ij += e.gpweight[igp] * e.vol * e.BaseGpVals[inode][igp] * e.BaseGpVals[jnode][igp];
        }
        c_ij *= h / k;
        C_coeffs.push_back( T( (*e.con)[inode], (*e.con)[jnode], c_ij ) );
      }
    }
  }
  if (isConvection) {
    lhsCoeffs.insert(lhsCoeffs.end(), C_coeffs.begin(), C_coeffs.end());
  }
}

void Problem::assembleConvectionRHS() {

  Eigen::VectorXd convectionRhs;
  convectionRhs.resize( domain.mesh->nnodes );
  convectionRhs.setZero();

  double k = material["k"];
  double h = material["h"];

  mesh::Element e;
  double c_i;

  // assemble
  vector<int> indicesConvectionFacets = convectionFacets.getTrueIndices();
  for ( int ifacet : indicesConvectionFacets ) {

    //TODO: Think about active elements!
    e = domain.getBoundaryFacet( ifacet );

    for (int inode = 0; inode < e.nnodes; ++inode) {
      c_i = 0.0;
      for (int igp = 0; igp < e.ngpoints; ++igp) {
        c_i += e.gpweight[igp] * e.vol * e.BaseGpVals[inode][igp] * Tenv;
      }
      c_i *= h / k;
      convectionRhs[(*e.con)[inode]] += c_i;
    }
  }
  rhs += convectionRhs;
}
