#include <iostream>
#include "Problem.h"

void Problem::assemble() {
  // ASSEMBLY
  // LHS & RHS, inactive nodes
  forceInactiveNodes();

  assembleSpatialPDE();// Spatial operator and source term

  assembleWeakBcs();

  // LHS & RHS, time
  assembleTime();

  // Proper assembly
  ls.assemble();

  // Dirichlet BC
  forceDirichletNodes();
}
