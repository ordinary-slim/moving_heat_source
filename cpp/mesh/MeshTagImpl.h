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
  x = vector<T>(  mesh->getNumEntities( dim ) , cte );
}

template<typename T>
MeshTag<T>::MeshTag(const mesh::Mesh *mesh, int dim, const std::vector<int> &indices){
  _mesh = mesh;
  _dim  = dim;
  x = vector<T>(  mesh->getNumEntities( dim ) , T(0) );
  for (int index : indices) {
    x[ index ] = T( 1 );
  }
}

template<typename T>
MeshTag<T>::MeshTag(const mesh::Mesh *mesh, std::vector<std::pair<int, T>> idxNvals, int dim){
  _mesh = mesh;
  _dim  = dim;
  x = vector<T>(  mesh->getNumEntities( dim ) , T(0) );
  for (std::pair<int, T> pair : idxNvals) {
    x[ pair.first ] = pair.second;
  }
}

template<typename T>
MeshTag<T>::MeshTag(const mesh::Mesh *mesh, const std::vector<int> &indices, const std::vector<T> &values, const int dim) {
  _mesh = mesh;
  _dim  = dim;
  x = vector<T>(  mesh->getNumEntities( dim ) , T(0) );
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
std::vector<int> MeshTag<T>::getIndices() const {
  /*
   * Get indices of MeshTag that evaluate to true.
   */
  return filterIndices( [](T tag){ return bool(tag);} );
}

template<typename T>
std::vector<int> MeshTag<T>::filterIndices( std::function<bool(T)> filter ) const {
  /*
   * Get indices of MeshTag that evaluate to filter[ MeshTag[index] ]
   */
  std::vector<int> indices;
  for (int index = 0; index < size(); ++index) {
    if ( filter( x[index] ) ) {
      indices.push_back( index );
    }
  }
  return indices;
}

template<typename T>
void MeshTag<T>::tag( std::function<bool(T)> filter, T value ) {
  /*
   * Tag indices which pass the filter
   */
  for (int index = 0; index < size(); ++index) {
    if (filter( x[index] ) ) {
      x[index] = value;
    }
  }
}
}
#endif
