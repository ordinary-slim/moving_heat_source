//https://stackoverflow.com/questions/10632251/undefined-reference-to-template-function
//The implementation of a non-specialized template must be visible to a translation unit that uses it.
#ifndef MESHTAGIMPL
#define MESHTAGIMPL
#include <stdexcept>
#include "MeshTag.h"
#include "Mesh.h"

namespace mesh {
template<typename T>
MeshTag<T>::MeshTag(const mesh::Mesh *mesh, int dim, const T &cte) {
  _mesh = mesh;
  _dim  = dim;
  x.resize( mesh->getNumEntities( dim ) );
  fill(x.begin(), x.end(), cte);
}

template<typename T>
MeshTag<T>::MeshTag(const mesh::Mesh *mesh, int dim, std::vector<T> values){
  _mesh = mesh;
  _dim  = dim;
  x.resize( mesh->getNumEntities( dim ) );
  if (values.size()==size()) {
    x = values;
  }
}

template<typename T>
MeshTag<T>::MeshTag(const mesh::Mesh *mesh, std::vector<std::pair<int, T>> idxNvals, int dim){
  _mesh = mesh;
  _dim  = dim;
  x.resize( mesh->getNumEntities( dim ) );
  for (std::pair<int, T> pair : idxNvals) {
    x[ pair.first ] = pair.second;
  }
}

template<typename T>
MeshTag<T>::MeshTag(const mesh::Mesh *mesh, const std::vector<int> &indices, const std::vector<T> &values, const int dim) {
  _mesh = mesh;
  _dim  = dim;
  x.resize( mesh->getNumEntities( dim ) );
  for (int subidx = 0; subidx < indices.size(); ++subidx) {
    x[ indices[subidx] ] = values[subidx];
  }
}

template<typename T>
T& MeshTag<T>::operator[](int idx) {
  return x[idx];
}

template<typename T>
const T& MeshTag<T>::operator[](int idx) const {
  return x[idx];
}

template<typename T>
void MeshTag<T>::setCteValue(const T &val) {
  // Set constant value
  fill(x.begin(), x.end(), val);
}
template<typename T>
void MeshTag<T>::setValues(std::vector<T> &v) {
  if (size()!=v.size() ) {
    throw std::invalid_argument( "Values provided don't match internal size." );
  } else {
    x = v;
  }
}

template<typename T>
std::vector<int> MeshTag<T>::getTrueIndices() {
  std::vector<int> indices;
  for (int index = 0; index < size(); ++index) {
    if ( bool( x[index] ) ) {
      indices.push_back( index );
    }
  }
  return indices;
}
}
#endif
