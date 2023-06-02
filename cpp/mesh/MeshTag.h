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

    MeshTag(const mesh::Mesh *mesh, int dim, const T &cte);

    MeshTag(const mesh::Mesh *mesh, int dim=0, std::vector<T> values= std::vector<T>());

    MeshTag(const mesh::Mesh *mesh, std::vector<std::pair<int, T>> idxNvals, int dim=0);

    MeshTag(const mesh::Mesh *mesh, const std::vector<int> &indices, const std::vector<T> &values, const int dim=0);

    T& operator[](int idx);

    const T& operator[](int idx) const;

    void setCteValue(const T &val);

    void setValues(std::vector<T> &v);

    std::vector<int> getIndices() const;
    std::vector<int> filterIndices( std::function<bool(T)> filter );

    int dim() const { return _dim; }
  private:
    int _dim;
    const Mesh* _mesh;
};
}
#endif
