#ifndef REFELEMENT
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "ElementTypes.h"
class ReferenceElement {
  public:
    int nnodes, dim, ngpoints;
    ElementType elementType;
    double vol = -1;
    Eigen::MatrixX3d pos, gpos;
    Eigen::Matrix3d XI_inverse;
    std::vector<std::vector<Eigen::Vector3d>> GradBaseGpVals;

    // Array of shape funcs
    std::vector<std::function<double(Eigen::Vector3d)>> shapeFuns;

    ReferenceElement(){}
    ReferenceElement( ElementType elType ) {
      elementType = elType;
      switch (elementType ) {
        case point1:
          // Not sure about this, necessary with structure of code
          nnodes = 1;
          ngpoints = 1;

          allocate( nnodes, ngpoints );

          dim =0;
          pos << 0.0, 0.0, 0.0;
          shapeFuns[0] = [](Eigen::Vector3d Xi) { return 1; };

          vol = 1;
          XI_inverse << 1.0, 0.0, 0.0,
                        0.0, 1.0, 0.0,
                        0.0, 0.0, 1.0;

          gpos = pos;//closed integration

          GradBaseGpVals[0][0] << 1.0, 0.0, 0.0;
          break;
        case line2:
          /*
           *  1 x___________x 2
          */
          nnodes = 2;
          ngpoints = 2;

          allocate( nnodes, ngpoints );

          dim =1;
          pos << -1.0, 0.0, 0.0,
                  1.0, 0.0, 0.0;
          shapeFuns[0] = [](Eigen::Vector3d Xi) { return 0.5*(1 - Xi(0) ); };
          shapeFuns[1] = [](Eigen::Vector3d Xi) { return 0.5*(1 + Xi(0) ); };

          vol = 2;
          XI_inverse << 0.5, 0.0, 0.0,
                        0.0, 1.0, 0.0,
                        0.0, 0.0, 1.0;

          gpos = pos;//closed integration

          for (int igp = 0; igp < ngpoints; igp++) {
            GradBaseGpVals[0][igp] << -0.5, 0.0, 0.0;
            GradBaseGpVals[1][igp] << +0.5, 0.0, 0.0;
          }
          break;
        case triangle3:
          /*
           * 3x
           *  |\
           *  | \
           *  |  \
           *  ^   \
           *  |    \
           * 1x__>__x2
           */
          nnodes = 3;
          ngpoints = 3;

          allocate( nnodes, ngpoints );

          dim =2;
          pos << 0.0, 0.0, 0.0,
                 1.0, 0.0, 0.0,
                 0.0, 1.0, 0.0;
          shapeFuns[0] = [](Eigen::Vector3d Xi) { return ( 1 - Xi(0) - Xi(1) ); };
          shapeFuns[1] = [](Eigen::Vector3d Xi) { return ( Xi(0) ); };
          shapeFuns[2] = [](Eigen::Vector3d Xi) { return ( Xi(1) ); };

          vol = 0.5;
          XI_inverse << +1.0, 0.0, 0.0,
                        0.0, +1.0, 0.0,
                        0.0, 0.0, +1.0;

          gpos = pos;//closed integration
                     //
          for (int igp = 0; igp < ngpoints; igp++) {
            GradBaseGpVals[0][igp] << -1.0, -1.0, 0.0;
            GradBaseGpVals[1][igp] << +1.0, 0.0, 0.0;
            GradBaseGpVals[2][igp] << 0.0, +1.0, 0.0;
          }
          break;
        case quad4:
          /*
           * 2x____________x1
           *  |            |
           *  |     ^      |
           *  |     |__>   |
           *  |            |
           *  |            |
           * 3x____________x4
           */
          nnodes = 4;
          ngpoints = 4;

          allocate( nnodes, ngpoints );

          dim =2;
          pos << 1.0, 1.0, 0.0,
                -1.0, 1.0, 0.0,
                -1.0,-1.0, 0.0,
                 1.0,-1.0, 0.0;
          shapeFuns[0] = [](Eigen::Vector3d Xi) {
            return ( 0.25*( 1 + Xi(0) )*( 1 + Xi(1)  ) );
          };
          shapeFuns[1] = [](Eigen::Vector3d Xi) {
            return ( 0.25*( 1 - Xi(0) )*( 1 + Xi(1) ) );
          };
          shapeFuns[2] = [](Eigen::Vector3d Xi) {
            return ( 0.25*( 1 - Xi(0) )*( 1 - Xi(1) ) );
          };
          shapeFuns[3] = [](Eigen::Vector3d Xi) {
            return ( 0.25*( 1 + Xi(0) )*( 1 - Xi(1) ) );
          };

          vol = 4;
          XI_inverse << -0.5, 0.5, 0.0,
                        0.0, -0.5, 0.0,
                        0.0, 0.0, 1.0;

          gpos = pos;//closed integration

          {
            //auto here is std::function<Eigen::Vector3d(Eigen::Vector3d)>
            auto GradBase1 = [](Eigen::Vector3d Xi) {
              return Eigen::Vector3d(+0.25*(1+Xi(1)), +0.25*(1+Xi(0)), 0.0);
            };
            auto GradBase2 = [](Eigen::Vector3d Xi) {
              return Eigen::Vector3d(-0.25*(1+Xi(1)), 0.25*(1-Xi(0)), 0.0);
            };
            auto GradBase3 = [](Eigen::Vector3d Xi) {
              return Eigen::Vector3d(-0.25*(1-Xi(1)), -0.25*(1-Xi(0)), 0.0);
            };
            auto GradBase4 = [](Eigen::Vector3d Xi) {
              return Eigen::Vector3d(0.25*(1-Xi(1)), -0.25*(1+Xi(0)), 0.0);
            };
            std::vector<std::function<Eigen::Vector3d(Eigen::Vector3d)>> GradBaseFuns = {GradBase1,
                                                                                        GradBase2,
                                                                                        GradBase3,
                                                                                        GradBase4};
            for (int inode = 0; inode < nnodes; inode++) {
              for (int igp = 0; igp < ngpoints; igp++) {
                GradBaseGpVals[inode][igp] = GradBaseFuns[inode]( gpos.row(igp) );
              }
            }
          break;
          }
          break;
        default:
          printf("Unknown element type\n");
          exit(EXIT_FAILURE);
      }
    }
  void allocate(int nnodes, int ngpoints ) {
    // allocate pos, shapeFuns, gpos, GradBaseGpVals
    pos.resize(nnodes, 3);
    shapeFuns.resize(nnodes);

    gpos.resize(ngpoints, 3);

    GradBaseGpVals.resize( nnodes );
    for (int inode = 0; inode < nnodes; inode++) {
      GradBaseGpVals[inode].resize( ngpoints );
    }
  }
};
#define REFELEMENT
#endif
