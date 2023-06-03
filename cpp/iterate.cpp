#include <iostream>
#include "Problem.h"

void Problem::iterate() {
  if (not(hasPreIterated)) {
    preIterate();
  }

  assemble();

  ls->solve();

  gather();

  postIterate();
}
