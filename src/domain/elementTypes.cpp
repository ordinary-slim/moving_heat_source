#include "elementTypes.h"
#include <iostream>

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
  std::cout << "Not implemented yet" << std::endl;
  exit(-1);
}

