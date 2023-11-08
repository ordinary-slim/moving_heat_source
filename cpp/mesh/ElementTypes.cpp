#include "ElementTypes.h"
#include <iostream>

static const char* nameofElType[] = {
  "point1", "line2", "triangle3", "quad4", "hexa8"
};

int getNnodesElType (ElementType elType ) {
  switch (elType) {
    case point1: return 1;
    case line2: return 2;
    case triangle3: return 3;
    case quad4: return 4;
    case hexa8: return 8;
    default: return -1;
  };
};

ElementType getFacetElType( ElementType Dent_elType ) {
  switch (Dent_elType) {
    case point1: return point1;
    case line2: return point1;
    case triangle3: return line2;
    case quad4: return line2;
    case hexa8: return quad4;
    default:
      throw std::invalid_argument("Not ready yet for this element type.");
  }
}

std::string elementType2String(ElementType elType) {
  return nameofElType[elType];
}

ElementType string2ElementType(std::string aux_cell_type) {
  //TODO: Move this to constructor
  //Is this possible for enum?
  ElementType cell_type_flag;
  if (aux_cell_type == "line2") {
    cell_type_flag = line2;
  } else if (aux_cell_type == "triangle3") {
    cell_type_flag = triangle3;
  } else if (aux_cell_type == "quad4") {
    cell_type_flag = quad4;
  } else if (aux_cell_type == "hexa8") {
    cell_type_flag = hexa8;
  } else {
    throw std::invalid_argument( "Not ready yet for this cell type." );
  }
  return cell_type_flag;
}

std::vector<std::vector<unsigned int>> getFacetVertexSets( const std::vector<unsigned int> &localCon_D0, ElementType entD_elType) {
  /*
   * Returns global indices of facets
   * Some spaghetti coding going on here
   * Gotta figure this out
   */

  std::vector<std::vector<unsigned int>> setsOfVertices;

  ElementType facetElType = getFacetElType( entD_elType );
  int nnodesCell = localCon_D0.size();
  int nnodesFacet = getNnodesElType(facetElType);

  std::vector<unsigned int> vset;//single set of vertices

  switch (entD_elType){
    case line2:
      setsOfVertices.push_back( std::vector<unsigned int>( 1, localCon_D0[0] ) );
      setsOfVertices.push_back( std::vector<unsigned int>( 1, localCon_D0[1] ) );
      break;
    case triangle3: case quad4:
      vset.resize( nnodesFacet );
      for (int ipoin = 0; ipoin < nnodesCell; ipoin++) {
        for (int jpoin = 0; jpoin < nnodesCell; jpoin++) {
          vset[jpoin] = localCon_D0[ (ipoin + jpoin)%nnodesCell ];
        }
        //std::sort( vset.begin(), vset.end() );
        setsOfVertices.push_back( vset );
      }
      break;
    case hexa8: {
      std::vector<std::vector<unsigned int>> setsOfIndices;
      setsOfIndices.push_back( std::vector<unsigned int>{0, 1, 2, 3} );
      setsOfIndices.push_back( std::vector<unsigned int>{4, 5, 6, 7} );
      setsOfIndices.push_back( std::vector<unsigned int>{0, 1, 5, 4} );
      setsOfIndices.push_back( std::vector<unsigned int>{3, 2, 6, 7} );
      setsOfIndices.push_back( std::vector<unsigned int>{0, 3, 7, 4} );
      setsOfIndices.push_back( std::vector<unsigned int>{1, 2, 6, 5} );
      for (auto& iset : setsOfIndices) {
        vset.clear();
        vset.reserve( nnodesCell );
        for (int index : iset) {
          vset.push_back( localCon_D0[ index ] );
        }
        setsOfVertices.push_back( vset );
      }
      break;
    }
    default:
      throw std::invalid_argument( "getFacetVertexSets not implemented yet for " + elementType2String(entD_elType) );
  }
  return setsOfVertices;
}

std::vector<std::vector<unsigned int>> getFacetVertexSets( ElementType entD_elType) {
  /*
   * Returns local indices of facets
   */
  std::vector<std::vector<unsigned int>> setsOfIndices;
  ElementType facetElType = getFacetElType( entD_elType );
  int nnodesCell  = getNnodesElType(entD_elType);
  int nnodesFacet = getNnodesElType(facetElType);
  std::vector<unsigned int> vset;//single set of vertices

  switch (entD_elType){
    case line2:
      setsOfIndices.push_back( std::vector<unsigned int>{0} );
      setsOfIndices.push_back( std::vector<unsigned int>{1} );
      break;
    case triangle3: case quad4:
      vset.resize( nnodesFacet );
      for (int ipoin = 0; ipoin < nnodesCell; ipoin++) {
        for (int jpoin = 0; jpoin < nnodesCell; jpoin++) {
          vset[jpoin] = (ipoin + jpoin)%nnodesCell;
        }
        setsOfIndices.push_back( vset );
      }
      break;
    case hexa8: {
      setsOfIndices.push_back( std::vector<unsigned int>{0, 1, 2, 3} );
      setsOfIndices.push_back( std::vector<unsigned int>{4, 5, 6, 7} );
      setsOfIndices.push_back( std::vector<unsigned int>{0, 1, 5, 4} );
      setsOfIndices.push_back( std::vector<unsigned int>{3, 2, 6, 7} );
      setsOfIndices.push_back( std::vector<unsigned int>{0, 3, 7, 4} );
      setsOfIndices.push_back( std::vector<unsigned int>{1, 2, 6, 5} );
      break;
    }
    default:
      throw std::invalid_argument( "getFacetVertexSets not implemented yet for " + elementType2String(entD_elType) );
  }

  return setsOfIndices;
}
