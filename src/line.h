#ifndef LINE
#include <iostream>
#include <vector>
#include "element.h"
using namespace std;

class Line: public Element {
  public:
    Line() {
      dimension = 1;
      nnodes = 2;
      rpos = {-1.0, +1.0};
      pos.reserve( 2 );
      // 2 nodes closed integration
      rgpos = {-1.0, +1.0};
      gpos.reserve( 2 );
    }

    void setBaseFunGradGpVals() {
      for (int jgp=0; jgp<nnodes; jgp++) {
        baseFunGradGpVals[0][jgp][0] = -1.0 / vol ;
        baseFunGradGpVals[1][jgp][0] = +1.0 / vol ;
      }
    }
};
#define LINE
#endif
