#include <iostream>
#include "Problem.h"

void Problem::iterate() {
  preIterate();

  assemble();

  ls->solve();

  gather();

  postIterate();
}
