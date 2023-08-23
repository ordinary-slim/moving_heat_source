#include <iostream>
#include "Problem.h"

void Problem::assemble(const Problem* externalProblem) {
  // ASSEMBLY
  assembleSpatialPDE();// Spatial operator and source term

  assembleWeakBcs();

  // LHS & RHS, time
  assembleTime();

  if (isCoupled and not( externalProblem == nullptr )) {
    assembleNeumannGamma( externalProblem );
    assembleDirichletGamma( externalProblem );
  }

  // If I own the LS, assemble it
  if (not(assembling2external)) {
    ls->assemble();
  }
}
