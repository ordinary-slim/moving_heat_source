#include "ElementTypes.h"
#include <iostream>

static const char* nameofElType[] = {"point1", "line2", "triangle3", "quad4"};

int getNnodesElType (ElementType elType ) {
  switch (elType) {
    case point1: return 1;
    case line2: return 2;
    case triangle3: return 3;
    case quad4: return 4;
    default: return -1;
  };
};

ElementType getIncidentElType( ElementType Dent_elType, int d ) {
  switch (d){
    case 0: return point1;
    case 1: return line2;
  }
  printf("%s not implemented yet for ", __FUNCTION__);
  std::cout << "( " << nameofElType[Dent_elType] << ", " << d << ")\n";
  exit(-1);
}

