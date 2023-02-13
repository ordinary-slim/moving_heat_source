#include "elementTypes.h"
#include <map>
#include <iostream>

std::map<ElementType, int> getNnodesElType {
  {point1, 1},
  {line2, 2},
  {triangle3, 3},
  {quad4, 4},
};

ElementType getDdElType( ElementType Dent_elType, int d ) {
  switch (d){
    case 0:
      return point1;
    case 1:
      return line2;
  }
  std::cout << "Not implemented yet" << std::endl;
  exit(-1);
}

