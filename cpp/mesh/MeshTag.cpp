#include "MeshTag.h"

namespace mesh {
  MeshTag<int> mark( const Mesh &mesh, int dim, const vector<int> &indices ) {
    vector<int> values = vector<int>( mesh.getNumEntities( dim ), 0 );
    for (int index : indices) {
      values[index] = 1;
    }
    return MeshTag<int>( &mesh, dim, values );
  }
}
