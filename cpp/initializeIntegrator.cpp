#include "problem.h"
#include <Eigen/Core>
#include "../external/pybind11/include/pybind11/eigen.h"

void Problem::initializeIntegrator(Eigen::MatrixXd pSols) {
  if (timeIntegrator.nstepsRequired > pSols.cols() ) {
    cout << "Not enough value provided for time integrator inititialization " << endl;
    exit(1);
  }
  unknown.prevValues = pSols(Eigen::placeholders::all, Eigen::seq( 0, timeIntegrator.nstepsRequired - 1));
  unknown.values = unknown.prevValues(Eigen::placeholders::all, 0);
  timeIntegrator.nstepsStored += pSols.cols();
}
