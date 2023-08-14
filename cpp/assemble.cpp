#include <iostream>
#include "Problem.h"

void Problem::assemble() {
  // ASSEMBLY
  assembleSpatialPDE();// Spatial operator and source term

  assembleWeakBcs();

  // LHS & RHS, time
  assembleTime();

  // If I own the LS, assemble it
  if (not(assembling2external)) {
    ls->assemble();
  }
}
