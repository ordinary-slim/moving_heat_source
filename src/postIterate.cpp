#include "includes/problem.h"
#include <Eigen/Core>

void Problem::postIterate() {
  // End iteration operations
  // Overwrite last column of prevSolutions
  prevSolutions.col( prevSolutions.cols()-1 ) << solution;
  // Permutate N-1, 0, ..., N-2
  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(timeIntegrator.nstepsRequired);
  Eigen::VectorXi indices(timeIntegrator.nstepsRequired);
  for (int i = 0; i < indices.size(); i++) {
    indices[i] = i-1;
  }
  indices[0] = indices.size() -  1 ;
  perm.indices() = indices;
  prevSolutions = prevSolutions * perm;

  ++timeIntegrator.nstepsStored;
  setTime( time + dt );
  ++iter;
}
