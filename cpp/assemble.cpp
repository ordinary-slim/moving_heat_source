#include <iostream>
#include "Problem.h"

void Problem::assemble() {
  if (not(assembling2external)) {
    ls->allocate();
  }
  // ASSEMBLY
  updateForcedDofs();

  assembleSpatialPDE();// Spatial operator and source term

  assembleWeakBcs();

  // LHS & RHS, time
  assembleTime();

  // Dirichlet BC
  forceDofs();

  // Proper assembly
  if (not(assembling2external)) {
    ls->assemble();
  }
}
