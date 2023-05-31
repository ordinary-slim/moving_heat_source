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
    //LHS
    vector<Eigen::Triplet<double>> timeDerivCoeffs;
    timeDerivCoeffs.resize( domain.massCoeffs.size() );
    for (int iMassEntry = 0; iMassEntry < domain.massCoeffs.size(); ++iMassEntry) {
      timeDerivCoeffs[iMassEntry] = Eigen::Triplet<double>(
          dofNumbering[domain.massCoeffs[iMassEntry].row()],
          dofNumbering[domain.massCoeffs[iMassEntry].col()],
          timeIntegrator.lhsCoeff*domain.massCoeffs[iMassEntry].value()/dt );
    }
    ls->lhsCoeffs.insert( ls->lhsCoeffs.end(), timeDerivCoeffs.begin(), timeDerivCoeffs.end() );

    //RHS
    int prevValCounter = 0;
    for (fem::Function prevFun: previousValues) {
      Eigen::VectorXd rhsContrib = domain.massMat * (prevFun.values * timeIntegrator.rhsCoeff[prevValCounter] ) / dt;
      for (int inode = 0; inode < domain.mesh->nnodes; ++inode) {
        ls->rhs[dofNumbering[inode]] += rhsContrib(inode);
      }
      ++prevValCounter;
    }
  }
}
