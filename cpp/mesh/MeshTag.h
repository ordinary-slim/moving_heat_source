#ifndef MESHTAG
#include <vector>
#include "Mesh.h"

namespace mesh {

template<typename T>
class MeshTag {
  public:
    std::vector<T> x;//values
    int size() const { return x.size(); }

    MeshTag(const mesh::Mesh *mesh, int dim=0, vector<T> values= vector<T>()){
      _mesh = mesh;
      _dim  = dim;
      x.resize( mesh->getNumEntities( dim ) );
      if (values.size()==size()) {
        x = values;
      }
    }

    MeshTag(const mesh::Mesh *mesh, vector<pair<int, T>> idxNvals, int dim=0){
      _mesh = mesh;
      _dim  = dim;
      x.resize( mesh->getNumEntities( dim ) );
      for (pair<int, T> pair : idxNvals) {
        x[ pair.first ] = pair.second;
      }
    }

    MeshTag(const mesh::Mesh *mesh, const vector<int> &indices, const vector<T> &values, const int dim=0){
      _mesh = mesh;
      _dim  = dim;
      x.resize( mesh->getNumEntities( dim ) );
      for (int subidx = 0; subidx < indices.size(); ++subidx) {
        x[ indices[subidx] ] = values[subidx];
      }
    }

    T& operator[](int idx) {
      return x[idx];
    }

    vector<int> getTrueIndices() {
      vector<int> indices;
      for (int index = 0; index < size(); ++index) {
        if ( bool( x[index] ) ) {
          indices.push_back( index );
        }
      }
      return indices;
    }

    void setValues(std::vector<T> &v) {
      if (size()!=v.size() ) {
        cout << "Bad call to MeshTag::setValues!" << endl;
        exit(-1);
      }
      x = v;
    }

    int dim() const { return _dim; }
  private:
    int _dim;
    const Mesh* _mesh;
};
MeshTag<int> mark( const Mesh &mesh, int dim = 0, const vector<int> &indices = vector<int>() );
}

#define MESHTAG
#endif
