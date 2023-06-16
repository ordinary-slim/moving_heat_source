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
    std::vector<unsigned int> refNodesMapping;
    Eigen::Matrix3d XI_inverse;
    std::vector<double> gpweight;
    std::vector<std::vector<double>> BaseGpVals;
    std::vector<std::vector<Eigen::Vector3d>> GradBaseGpVals;
    bool openIntegration = false;//default closed integration

    // Array of shape funcs
    std::vector<std::function<double(Eigen::Vector3d)>> shapeFuns;
    std::vector<std::function<Eigen::Vector3d(Eigen::Vector3d)>> gradShapeFuns;

    ReferenceElement(){}
    ReferenceElement( ElementType elType, int ngps = -1 );
    void allocate(int nnodes, int ngpoints );
};
#endif
