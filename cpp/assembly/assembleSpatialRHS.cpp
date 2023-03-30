#include "../Problem.h"
#include <iostream>
#include <Eigen/Core>
#include <cmath>
#include <vector>

void Problem::assembleSpatialRHS() {
  mhs.pulse.resize( mesh.nnodes ); // source term
  mhs.pulse.setZero();

  double r_i;
  Eigen::Vector3d x_gp;

  mesh::Element e;
  // assemble
  for (int ielem = 0; ielem < mesh.nels; ++ielem) {
    if (mesh.activeElements[ielem]==0){
      continue;
    }
    e = mesh.getElement( ielem );
    for (int inode = 0; inode < e.nnodes; ++inode) {
      r_i = 0;
      for (int igp = 0; igp < e.ngpoints; ++igp) {
        x_gp = e.gpos.row( igp );
        r_i += e.gpweight[igp] * e.BaseGpVals[inode][igp] * e.vol *
          mhs.powerDensity(x_gp, time, mhs.currentPosition, mhs.power, mhs.efficiency, mhs.radius);
      }
      mhs.pulse[e.con[inode]] += r_i;
    }
  }
  rhs += mhs.pulse;

  // DEBUGGING: For sending pulse to post
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  //Solve linear system
  solver.compute( M );
  if (not(solver.info() == Eigen::Success)) {
    std::cout << "Singular matrix!" << std::endl;
  }
  mhs.pulse = solver.solve(mhs.pulse);
}
