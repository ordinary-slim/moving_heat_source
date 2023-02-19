#include "problem.h"
#include <Eigen/Core>

void Problem::postIterate() {
  // End iteration operations
  // Overwrite last column of unknown.prevValues
  unknown.prevValues.col( unknown.prevValues.cols()-1 ) << unknown.values;
  // Permutate N-1, 0, ..., N-2
  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(timeIntegrator.nstepsRequired);
  Eigen::VectorXi indices(timeIntegrator.nstepsRequired);
  for (int i = 0; i < indices.size(); i++) {
    indices[i] = i-1;
  }
  indices[0] = indices.size() -  1 ;
  perm.indices() = indices;
  unknown.prevValues = unknown.prevValues * perm;

  ++timeIntegrator.nstepsStored;
  setTime( time + dt );
  ++iter;
}
