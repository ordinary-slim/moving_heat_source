#ifndef CONNEC
#include <Eigen/Core>
#include "ElementTypes.h"
#include "RefElement.h"

namespace mesh
{
class Connectivity {
  /*
   * Connectivity mesh entitites oDim -> tDim
  */
  public:
    int oDim = -1;
    int tDim = -1;
    int nels_oDim = -1;
    int nels_tDim = -1;
    ElementType oelType;
    ElementType telType;
    std::vector<std::vector<unsigned int>> con;// connectivy

    Connectivity() {}

    Connectivity( std::vector<std::vector<unsigned int>> in_con, int in_oDim, int in_tDim, int in_nels_oDim, int in_nels_tDim,
        ElementType in_oelType, ElementType in_telType) {
      con = in_con;
      oDim = in_oDim;
      tDim = in_tDim;
      nels_oDim = in_nels_oDim;
      nels_tDim = in_nels_tDim;
      oelType = in_oelType;
      telType = in_telType;
    }

    const std::vector<unsigned int>* getLocalCon( int idx ) const {
      return &con[ idx ];
    }
    /*
    std::vector<int> getLocalConVector( int idx ) {
      return std::vector<int>( con.row(idx).data(), con.row(idx).data() + con.cols() );
    }
    */
};

Connectivity transpose(Connectivity inCon);
Connectivity intersect(Connectivity inCon1, Connectivity inCon2);
std::tuple<Connectivity, Connectivity> buildBoundaryConnectivities(
    Connectivity DO_connec, Connectivity DD_connec);

}
#define CONNEC
#endif
