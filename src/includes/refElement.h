#ifndef REFELEMENT
#include <iostream>
#include <vector>
#include <Eigen/Dense>
class refElement {
  public:
    int nnodes, dim, ngpoints;
    int elementType;
    double vol = -1;
    Eigen::MatrixX3d pos;
    Eigen::Matrix3d XI_inverse;
    std::vector<std::vector<Eigen::Vector3d>> GradBaseGpVals;

    refElement(){}
    refElement( int elType ) {
      elementType = elType;
      switch (elementType ) {
        case 0://line2
          /*
           *  1 x___________x 2
          */
          nnodes = 2;
          ngpoints = 2;
          dim =1;
          pos.resize(nnodes, 3);
          pos << -1.0, 0.0, 0.0,
                  1.0, 0.0, 0.0;
          vol = 2;
          XI_inverse << 0.5, 0.0, 0.0,
                        0.0, 1.0, 0.0,
                        0.0, 0.0, 1.0;

          GradBaseGpVals.resize( nnodes );
          GradBaseGpVals[0].resize( ngpoints );
          GradBaseGpVals[1].resize( ngpoints );
          for (int igp = 0; igp < ngpoints; igp++) {
            GradBaseGpVals[0][igp] << -0.5, 0.0, 0.0;
            GradBaseGpVals[1][igp] << +0.5, 0.0, 0.0;
          }
          break;
        case 3://triangle3
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
          dim =2;
          pos.resize(nnodes, 3);
          pos << 0.0, 0.0, 0.0,
                 1.0, 0.0, 0.0,
                 0.0, 1.0, 0.0;
          vol = 1;
          XI_inverse << +1.0, 0.0, 0.0,
                        0.0, +1.0, 0.0,
                        0.0, 0.0, +1.0;
          break;
        case 4://quad4
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
          dim =2;
          pos.resize(nnodes, 3);
          pos << 1.0, 1.0, 0.0,
                -1.0, 1.0, 0.0,
                -1.0,-1.0, 0.0,
                 1.0,-1.0, 0.0;
          vol = 4;
          XI_inverse << -0.5, 0.5, 0.0,
                        0.0, -0.5, 0.0,
                        0.0, 0.0, 1.0;
          break;
        default:
          printf("Unknown element type\n");
          exit(EXIT_FAILURE);
      }
    }
};
#define REFELEMENT
#endif
