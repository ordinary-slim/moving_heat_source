#include <iostream>
#include "Problem.h"

void Problem::assemble() {
  if (not(assembling2external)) {
    ls->allocate();
  }
  // ASSEMBLY
  // LHS & RHS, inactive nodes
  forceInactiveNodes();

  assembleSpatialPDE();// Spatial operator and source term

  assembleWeakBcs();

  // LHS & RHS, time
  assembleTime();

  // Dirichlet BC
  updateForcedDofs();

  // Proper assembly
  if (not(assembling2external)) {
    ls->assemble();
  }
}
