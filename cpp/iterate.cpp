#include <iostream>
#include "Problem.h"

void Problem::iterate() {
  if (not(hasPreIterated)) {
    preIterate(true);
  }

  assemble();

  ls->solve();

  gather();

  postIterate();
}
