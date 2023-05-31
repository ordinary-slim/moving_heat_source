#include <iostream>
#include "Problem.h"

void Problem::iterate() {
  preIterate();

  assemble();

  ls.solve();
  unknown.values = ls.sol;

  postIterate();
}
