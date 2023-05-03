#ifndef ELTYPE
#include <string>
enum ElementType {
  point1,
  line2,
  triangle3,
  quad4,
};
int getNnodesElType( ElementType );
ElementType getIncidentElType( ElementType Dent_elType, int d );
std::string getElementTypeName(ElementType elType);
#define ELTYPE
#endif
