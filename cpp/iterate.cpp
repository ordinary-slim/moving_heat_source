#include <iostream>
#include "Problem.h"

void Problem::iterate() {
  preIterate();

  assemble();

  solve();

  postIterate();
}
