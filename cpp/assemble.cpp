#include <iostream>
#include "Problem.h"

void Problem::assemble() {
  // ASSEMBLY
  // LHS & RHS, inactive nodes
  forceInactiveNodes();
  // LHS, space
  assembleSpatialPDE();

  // Neumann BC
  assembleNeumann();

  // Convection BC
  assembleConvectionLHS();
  assembleConvectionRHS();

  // LHS & RHS, VMS
  //if (isStabilized) {
    //assembleStabilization();
  //}

  // LHS & RHS, time
  assembleTime();

  // Proper assembly
  lhs.setFromTriplets( lhsCoeffs.begin(), lhsCoeffs.end() );

  // Dirichlet BC
  forceDirichletNodes();
}
