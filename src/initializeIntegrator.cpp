#include "problem.h"
#include <Eigen/Core>
#include "../external/pybind11/include/pybind11/eigen.h"

void Problem::initializeIntegrator(Eigen::MatrixXd pSols) {
  cout << "prevSolutions, before = " << prevSolutions << endl;
  prevSolutions = pSols;
  cout << "pSols = " << pSols << endl;
  cout << "prevSolutions, after = " << prevSolutions << endl;
  nstepsStored += pSols.cols();
}
