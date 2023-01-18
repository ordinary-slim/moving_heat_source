#include "../includes/mesh.h"


void Mesh::initialize1DMesh( double A, double B, int numberOfEls ){
  nels = numberOfEls;
  nnodes = nels + 1;
  double L = (B - A);
  double h = L / nels;
  cout << "h = " << h << endl;

  con.resize( nels );
  elementTypes.resize( nels );
  pos.resize( nnodes );

  for (int i = 0; i < nnodes; i++) {
    pos[i] = A + i*h;
  }

  // build connectivity; inneficient but clear
  nnodes_per_el = 2;
  for (int i = 0; i < nels; i++) {
    con[i].push_back( i );
    con[i].push_back( i+1 );
    elementTypes[i] = 0;
  }
}
