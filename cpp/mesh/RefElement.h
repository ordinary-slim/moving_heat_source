#ifndef REFELEMENT
#define REFELEMENT
#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include "ElementTypes.h"
class ReferenceElement {
  public:
    int nnodes, dim, ngpoints = -1;
    ElementType elementType;
    double vol = -1;
    Eigen::MatrixX3d pos, gpos;
    Eigen::Matrix3d XI_inverse;
    std::vector<double> gpweight;
    std::vector<std::vector<double>> BaseGpVals;
    std::vector<std::vector<Eigen::Vector3d>> GradBaseGpVals;
    bool openIntegration = false;//default closed integration

    // Array of shape funcs
    std::vector<std::function<double(Eigen::Vector3d)>> shapeFuns;
    std::vector<std::function<Eigen::Vector3d(Eigen::Vector3d)>> gradShapeFuns;

    ReferenceElement(){}
    ReferenceElement( ElementType elType, int ngps = -1 ) {
      elementType = elType;
      nnodes = getNnodesElType( elType );
      if (ngps >= nnodes ) {
        ngpoints = ngps;
        openIntegration = true;
      } else {
        ngpoints = nnodes;//closed integration
      }
      allocate( nnodes, ngpoints );

      switch (elementType ) {
        case point1:
          // Not sure about this, necessary with structure of code

          dim =0;
          pos << 0.0, 0.0, 0.0;
          shapeFuns[0] = [](Eigen::Vector3d Xi) { return 1; };

          vol = 1;
          XI_inverse << 1.0, 0.0, 0.0,
                        0.0, 1.0, 0.0,
                        0.0, 0.0, 1.0;

          if (openIntegration) {
            switch (ngpoints) {
              default:
                printf("Open integration with %i not implemented for point1 element\n", ngpoints);
                exit(-1);
            }
          } else {
            gpos = pos;//closed integration
            fill( gpweight.begin(), gpweight.end(), 1.0 / nnodes );
          }

          gradShapeFuns[0] = [](Eigen::Vector3d Xi) { return Eigen::Vector3d(+1.0, 0.0, 0.0); };
          break;
        case line2:
          /*
           *  1 x___________x 2
          */

          dim =1;
          pos << -1.0, 0.0, 0.0,
                  1.0, 0.0, 0.0;
          shapeFuns[0] = [](Eigen::Vector3d Xi) { return 0.5*(1 - Xi(0) ); };
          shapeFuns[1] = [](Eigen::Vector3d Xi) { return 0.5*(1 + Xi(0) ); };

          vol = 2;
          XI_inverse << 0.5, 0.0, 0.0,
                        0.0, 1.0, 0.0,
                        0.0, 0.0, 1.0;

          if (openIntegration) {
            switch (ngpoints) {
              case 2:
                gpos << -1.0 / sqrt(3.0), 0.0, 0.0,
                        +1.0 / sqrt(3.0), 0.0, 0.0;
                gpweight = {0.5, 0.5};
                break;
              case 3:
                gpos << 0.0, 0.0, 0.0,
                        -sqrt( 3.0 / 5.0 ), 0.0, 0.0,
                        +sqrt( 3.0 / 5.0 ), 0.0, 0.0;
                gpweight = {8.0/18.0, 5.0/18.0, 5.0/18.0};
                break;
              default:
                printf("Open integration with %i nodes not implemented for line2 element.\n", ngpoints);
                exit(-1);
            }
          } else {
            gpos = pos;//closed integration
            fill( gpweight.begin(), gpweight.end(), 1.0 / nnodes );
          }

          gradShapeFuns[0] = [](Eigen::Vector3d Xi) { return Eigen::Vector3d(-0.5, 0.0, 0.0); };
          gradShapeFuns[1] = [](Eigen::Vector3d Xi) { return Eigen::Vector3d(+0.5, 0.0, 0.0); };
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

          if (openIntegration) {
            switch (ngpoints) {
              default:
                printf("Open integration with %i not implemented for this element\n", ngpoints);
                exit(-1);
            }
          } else {
            gpos = pos;//closed integration
            fill( gpweight.begin(), gpweight.end(), 1.0 / nnodes );
          }

          gradShapeFuns[0] = [](Eigen::Vector3d Xi) { return Eigen::Vector3d(-1.0, -1.0, 0.0); };
          gradShapeFuns[1] = [](Eigen::Vector3d Xi) { return Eigen::Vector3d(+1.0, 0.0, 0.0); };
          gradShapeFuns[2] = [](Eigen::Vector3d Xi) { return Eigen::Vector3d(0.0, +1.0, 0.0); };

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

          if (openIntegration) {
            switch (ngpoints) {
              case 4: {
                double invSqrt3 = 1 / sqrt( 3.0 );
                gpos << -invSqrt3, -invSqrt3, 0.0,
                        -invSqrt3, +invSqrt3, 0.0,
                        +invSqrt3, +invSqrt3, 0.0,
                        +invSqrt3, -invSqrt3, 0.0;
                gpweight = {1.0/4.0, 1.0/4.0, 1.0/4.0, 1.0/4.0};
                break;}
              default:
                printf("Open integration with %i not implemented for this element\n", ngpoints);
                exit(-1);
            }
          } else {
            gpos = pos;//closed integration
            fill( gpweight.begin(), gpweight.end(), 1.0 / nnodes );
          }

          gradShapeFuns[0] = [](Eigen::Vector3d Xi) {
            return Eigen::Vector3d(+0.25*(1+Xi(1)), +0.25*(1+Xi(0)), 0.0);
          };
          gradShapeFuns[1] = [](Eigen::Vector3d Xi) {
            return Eigen::Vector3d(-0.25*(1+Xi(1)), 0.25*(1-Xi(0)), 0.0);
          };
          gradShapeFuns[2] = [](Eigen::Vector3d Xi) {
            return Eigen::Vector3d(-0.25*(1-Xi(1)), -0.25*(1-Xi(0)), 0.0);
          };
          gradShapeFuns[3] = [](Eigen::Vector3d Xi) {
            return Eigen::Vector3d(0.25*(1-Xi(1)), -0.25*(1+Xi(0)), 0.0);
          };
          break;
        default:
          printf("Unknown element type\n");
          exit(EXIT_FAILURE);
      }

      // Compute BaseGpVals
      for (int inode = 0; inode < nnodes; ++inode){
        for (int igpoin = 0; igpoin < ngpoints; ++igpoin){
          BaseGpVals[inode][igpoin] = shapeFuns[inode]( gpos.row( igpoin ) );
        }
      }
      // Compute GradBaseGpVals
      for (int inode = 0; inode < nnodes; inode++) {
        for (int igp = 0; igp < ngpoints; igp++) {
          GradBaseGpVals[inode][igp] = gradShapeFuns[inode]( gpos.row(igp) );
        }
      }
    }
  void allocate(int nnodes, int ngpoints ) {
    // allocate pos, shapeFuns, gpos, GradBaseGpVals
    pos.resize(nnodes, 3);
    shapeFuns.resize(nnodes);
    gradShapeFuns.resize(nnodes);

    gpos.resize(ngpoints, 3);
    // Quadrature weights
    gpweight.resize( ngpoints );
    // BaseFun
    BaseGpVals.resize( nnodes );
    for (int inode = 0; inode < nnodes; ++inode) {
      BaseGpVals[inode].resize( ngpoints );
    }
    // GradBaseFun
    GradBaseGpVals.resize( nnodes );
    for (int inode = 0; inode < nnodes; inode++) {
      GradBaseGpVals[inode].resize( ngpoints );
    }
  }
};
#endif
