#include <iostream>
#include "Problem.h"
#include "linearAlgebra/LinearSystem.h"

void Problem::assemble(const Problem* externalProblem) {
  // ASSEMBLY
  assembleDomain();// Spatial operator and source term

  assembleBoundary();

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
