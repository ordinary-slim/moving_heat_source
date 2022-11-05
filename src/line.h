#ifndef LINE
#include <iostream>
#include <vector>
#include "element.h"
using namespace std;

class Line: public Element {
  public:
    Line() {
      nnodes = 2;
      rpos = {-1.0, +1.0};
      pos.reserve( 2 );
      // 2 nodes closed integration
      rgpos = {-1.0, +1.0};
      gpos.reserve( 2 );
    }
};
#define LINE
#endif
