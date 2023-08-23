#ifndef MESHTAG
#define MESHTAG
#include <vector>
#include <utility>
#include <functional>

namespace mesh {

class Mesh;

template<typename T>
class MeshTag {
  public:
    std::vector<T> x;//values
    int size() const { return x.size(); }

    MeshTag(const mesh::Mesh *mesh, int dim = 0, const T &cte = 0);

    MeshTag(const mesh::Mesh *mesh, int dim, const std::vector<int> &indices);

    MeshTag(const mesh::Mesh *mesh, std::vector<std::pair<int, T>> idxNvals, int dim=0);

    MeshTag(const mesh::Mesh *mesh, const std::vector<int> &indices, const std::vector<T> &values, const int dim=0);

    MeshTag(const mesh::MeshTag<T> &tag) = default;

    T& operator[](int idx);

    const T& operator[](int idx) const;

    void setCteValue(const T &val);

    std::vector<int> getIndices() const;
    std::vector<int> filterIndices( std::function<bool(T)> filter ) const;
    void tag( std::function<bool(T)> filter, T value );

    MeshTag<int>& operator&=( const MeshTag<T>& rhs ) {
      if ( (_mesh != rhs._mesh) || (_dim != rhs._dim ) ) {
        throw std::invalid_argument("Incompatible MeshTags");
      }
      for (int ient = 0; ient < size(); ++ient) {
        x[ient] = int( x[ient] & rhs[ient] );
      }
      return *this;
    }


    friend MeshTag<int> operator!(const MeshTag<T> &tag) {
      MeshTag<T> complement = MeshTag<T>( tag );//copy input tag
      for (int ient = 0; ient < complement.size(); ++ient) {
        complement[ient] = int( not( complement[ient] ) );
      }
      return complement;
    }

    int dim() const { return _dim; }

  private:
    int _dim;
    const Mesh* _mesh;
};
}
#endif
