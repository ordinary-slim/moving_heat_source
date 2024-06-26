#include "../linearAlgebra/LinearSystem.h"
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
    vector<Eigen::Triplet<double>> discreteTimeDerivCoeffs;
    discreteTimeDerivCoeffs.reserve( timeDerivCoeffs.size() );
    for (int indexContrib = 0; indexContrib < timeDerivCoeffs.size(); ++indexContrib) {
      int inodeDof = freeDofsNumbering[ timeDerivCoeffs[indexContrib].row() ];
      if ( inodeDof < 0 ) { continue; }// if forced node, keep going
      int jnodeGlobal = timeDerivCoeffs[indexContrib].col();
      int jnodeDof = dofNumbering[ jnodeGlobal ];
      double coeff = timeIntegrator.lhsCoeff*timeDerivCoeffs[indexContrib].value()/dt;
      if ( jnodeDof < 0 ) {
        // To RHS
        ls->rhs[inodeDof] += - coeff * unknown.values[ jnodeGlobal ];
      } else {
        // To LHS
        discreteTimeDerivCoeffs.push_back(  Eigen::Triplet<double>(
            inodeDof, jnodeDof, coeff) );
      }
    }
    ls->lhsCoeffs.insert( ls->lhsCoeffs.end(), discreteTimeDerivCoeffs.begin(), discreteTimeDerivCoeffs.end() );

    //RHS
    int prevValCounter = 0;
    for (fem::Function& prevFun: previousValues) {
      Eigen::VectorXd rhsContrib = timeDerivMat *
          (prevFun.values * timeIntegrator.rhsCoeff[prevValCounter] ) / dt;
      for (int inode = 0; inode < domain.mesh->nnodes; ++inode) {
        int inodeDof = freeDofsNumbering[inode];
        if (inodeDof >= 0) {
          ls->rhs[inodeDof] += rhsContrib(inode);
        }
      }
      ++prevValCounter;
    }
  }
}
