#include "../Problem.h"

void Problem::assembleTime() {
  // general treatment implicit schemes
  if (not isSteady) {
    //Set time integration
    if (timeIntegrator.nstepsStored < timeIntegrator.nstepsRequired ) {
      if (timeIntegrator.nstepsStored >= 4) {
        timeIntegrator.setCurrentIntegrator( 4 );
      } else if (timeIntegrator.nstepsStored >= 3) {
        timeIntegrator.setCurrentIntegrator( 3 );
      } else if (timeIntegrator.nstepsStored >= 2) {
        timeIntegrator.setCurrentIntegrator( 2 );
      } else {
        timeIntegrator.setCurrentIntegrator( 1 );
      }
    } else {
      timeIntegrator.setCurrentIntegrator( timeIntegrator.desiredIntegrator );
    }
    //Add time dependency
    vector<Eigen::Triplet<double>> timeDerivCoeffs;
    timeDerivCoeffs.resize( massCoeffs.size() );
    for (int iMassEntry = 0; iMassEntry < massCoeffs.size(); ++iMassEntry) {
      timeDerivCoeffs[iMassEntry] = Eigen::Triplet<double>( massCoeffs[iMassEntry].row(),
                                                            massCoeffs[iMassEntry].col(),
                                    timeIntegrator.lhsCoeff*massCoeffs[iMassEntry].value()/dt );
    }
    lhsCoeffs.insert( lhsCoeffs.end(), timeDerivCoeffs.begin(), timeDerivCoeffs.end() );

    int prevValCounter = 0;
    for (fem::Function prevFun: previousValues) {
      rhs += M * (prevFun.values * timeIntegrator.rhsCoeff[prevValCounter] ) / dt;
      ++prevValCounter;
    }
  }
}
