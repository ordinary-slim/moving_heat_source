#ifndef LUMPEDHS
#define LUMPEDHS
#include "HeatSource.h"
#include "../additiveManufacturing/Hatch.h"
#include "../mesh/ActiveMesh.h"

namespace heat {
class LumpedHeatSource : public HeatSource {
  private:
    HatchCollider hc;
  public:
    mesh::MeshTag<int> heatedElements;
    mesh::MeshTag<double> elementPulse;//DEBUG
    double pd;
    double heatedVolume;

    LumpedHeatSource( double width, double height,
        pybind11::dict &input, Problem *problem );

    void markHeatedElements( const Eigen::Vector3d &p1, const Eigen::Vector3d &p2 );
    void computeHeatedVolume();
    void computePowerDensity();

};
}
#endif
