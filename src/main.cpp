#include <iostream>
#include <vector>
#include "line.h"
#include "mesh.h"
using namespace std;

int main() {
  // generate mesh
  // initialize solution
  // loop
    // timestep
    // plot
    //
  // Buid domain + mesh
  double L, R, h;
  int nels;

  L = 1;
  nels = 10;

  Mesh mesh;
  mesh.initialize1DMesh( 0.0, L, nels);

  // Assembly
  Line l;
  for (int ielem = 0; ielem < mesh.nels; ielem++ ) {
    l = mesh.getElement( ielem );
    l.print();
  }
}
