#ifndef ELTYPE
#define ELTYPE
#include <vector>
#include <algorithm>
#include <string>
enum ElementType {
  point1,
  line2,
  triangle3,
  quad4,
  hexa8,
};
int getNnodesElType( ElementType );
ElementType getFacetElType( ElementType Dent_elType );
std::string elementType2String(ElementType elType);
ElementType string2ElementType(std::string aux_cell_type);
std::vector<std::vector<unsigned int>> getFacetVertexSets(
    const std::vector<unsigned int> &localCon_D0, ElementType entD_elType);
#endif
