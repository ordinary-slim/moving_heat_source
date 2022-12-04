#include "problem.h"
#include <Eigen/Core>
#include "../external/pybind11/include/pybind11/eigen.h"

void Problem::initializeIntegrator(Eigen::MatrixXd pSols) {
  if (nstepsRequired > pSols.cols() ) {
    cout << "Not enough value provided for time integrator inititialization " << endl;
    exit(1);
  }
  prevSolutions = pSols(Eigen::placeholders::all, Eigen::seq( 0, nstepsRequired - 1));
  cout << "prevSolutions=" << prevSolutions;
  nstepsStored += pSols.cols();
}
