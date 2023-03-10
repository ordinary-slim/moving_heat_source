#include "Problem.h"
#include "FEMFunction.h"
#include "mesh/Mesh.h"
#include "../external/pybind11/include/pybind11/eigen.h"
#include <Eigen/Core>

void Problem::updateFRFpos() {
  // Pre-iteration operations

  // update positions in no advection RF
  // done in pre iterate because activation is also done here
  mesh.shiftFRF += dt * mesh.speedFRF;
  for (int inode=0; inode < mesh.nnodes; inode++){
    mesh.posFRF.row( inode ) += dt * mesh.speedFRF;
  }
}

void Problem::preIterate() {
  /* Beginning of iteration operations*/
  // CLEANUP Linear System
  lhs.resize( mesh.nnodes, mesh.nnodes );
  rhs.resize( mesh.nnodes );
  lhs.setZero();
  rhs.setZero();
  // UPDATE to tn+1
  mhs.updatePosition( dt );
  setTime( time + dt );
  ++iter;
}


void Problem::postIterate() {
  /* End iteration operations */
  // STORE last timestep for time-integratino
  // Overwrite last column of unknown.prevValues
  unknown.prevValues.col( unknown.prevValues.cols()-1 ) << unknown.values;
  // Permutate N-1, 0, ..., N-2
  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(timeIntegrator.nstepsRequired);
  Eigen::VectorXi indices(timeIntegrator.nstepsRequired);
  indices[0] = indices.size() -  1 ;
  for (int i = 1; i < indices.size(); i++) {
    indices[i] = i-1;
  }
  perm.indices() = indices;
  unknown.prevValues = unknown.prevValues * perm;

  ++timeIntegrator.nstepsStored;
}

void Problem::initializeIntegrator(Eigen::MatrixXd pSols) {
  unknown.prevValues(Eigen::placeholders::all, Eigen::seq( 0, pSols.cols()-1)) = pSols;
  unknown.values = unknown.prevValues(Eigen::placeholders::all, 0);
  timeIntegrator.nstepsStored = pSols.cols();
}
