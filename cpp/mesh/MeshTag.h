#ifndef MESHTAG
#include <vector>
#include "Mesh.h"

namespace mesh {

template<typename T>
class MeshTag {
  public:
    std::vector<T> x;//values
    int size() const { return _size; }

    MeshTag() = default;

    MeshTag(mesh::Mesh *mesh, int dim=0, vector<T> values= vector<T>()){
      _mesh = mesh;
      _dim  = dim;
      if (_dim == _mesh->dim ) {
        _size = _mesh->nels;
      } else if ( _dim == _mesh->dim-1 ) {
        _size = _mesh->con_FacetCell.nels_oDim;
      } else if ( _dim == 0 ) {
        _size = _mesh->nnodes;
      } else {
        cout << "Not ready yet!" << endl;
        exit(-1);
      }
      x.resize( _size );
      if (values.size()==_size) {
        x = values;
      }
    }

    vector<int> getTrueIndices() {
      vector<int> indices;
      for (int index = 0; index < _size; ++index) {
        if ( bool( x[index] ) ) {
          indices.push_back( index );
        }
      }
      return indices;
    }

    void setValues(std::vector<T> &v) {
      x = v;
    }
    int dim() const { return _dim; }
  private:
    int _dim, _size;
    Mesh* _mesh;
};
}

#define MESHTAG
#endif
