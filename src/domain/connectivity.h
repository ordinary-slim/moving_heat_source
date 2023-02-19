#ifndef CONNEC
#include <Eigen/Core>
#include "elementTypes.h"
#include "refElement.h"

namespace mesh
{
class Connectivity {
  public:
    int oDim = -1;
    int tDim = -1;
    int nels_oDim = -1;
    int nels_tDim = -1;
    ElementType oelType;
    ElementType telType;
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic,
      Eigen::RowMajor> con;// connectivy

    Connectivity() {}

    Connectivity( Eigen::MatrixXi in_con, int in_oDim, int in_tDim, int in_nels_oDim, int in_nels_tDim,
        ElementType in_oelType, ElementType in_telType) {
      con = in_con;
      oDim = in_oDim;
      tDim = in_tDim;
      nels_oDim = in_nels_oDim;
      nels_tDim = in_nels_tDim;
      oelType = in_oelType;
      telType = in_telType;
    }

    Eigen::VectorXi getLocalCon( int idx ) {
      return con.row( idx );
    }
};

Connectivity transpose(Connectivity inCon);
Connectivity intersect(Connectivity inCon1, Connectivity inCon2);
std::tuple<Connectivity, Connectivity> build(int d, Connectivity DO_connec, Connectivity DD_connec);

}
#define CONNEC
#endif
