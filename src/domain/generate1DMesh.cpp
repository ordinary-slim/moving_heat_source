#include "mesh.h"
#include "element.h"


void mesh::Mesh::generate1DMesh( double A, double B, int numberOfEls ){
  nels = numberOfEls;
  nnodes = nels + 1;
  double L = (B - A);
  double h = L / nels;
  cout << "h = " << h << endl;

  con_CellPoint.con.resize( nels, 2 );
  elementTypes.resize( nels );
  pos.resize( nnodes, 3 );

  for (int i = 0; i < nnodes; i++) {
    pos(i, 0) = A + i*h;
    pos(i, 1) = 0.0;
    pos(i, 2) = 0.0;
  }

  // build connectivity; inneficient but clear
  nnodes_per_el = 2;
  for (int i = 0; i < nels; i++) {
    con_CellPoint.con(i, 0 ) = i;
    con_CellPoint.con(i, 1 ) = i+1;
    elementTypes[i] = line2;
  }
}
