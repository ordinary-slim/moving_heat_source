#ifndef LUMPEDHS
#define LUMPEDHS
#include "HeatSource.h"
#include "../additiveManufacturing/Hatch.h"
#include "../mesh/ActiveMesh.h"

namespace heat {
class LumpedHeatSource : public HeatSource {
  private:
    HatchCollider hc;
    mesh::ActiveMesh *domain;
  public:
    mesh::MeshTag<int> heatedElements;
    mesh::MeshTag<double> elementPulse;//DEBUG
    double pd;
    double heatedVolume;

    LumpedHeatSource( mesh::ActiveMesh *domain, double width, double height ) :
      hc( domain->mesh, width, height ),
      heatedElements( mesh::MeshTag<int>( domain->mesh, domain->mesh->dim, 0 ) ),
      elementPulse( mesh::MeshTag<double>( domain->mesh, domain->mesh->dim, 0.0 ) )//DEBUG
    {
      type = lumped;
      this->domain = domain;
    }

    void markHeatedElements( const Eigen::Vector3d &p1, const Eigen::Vector3d &p2 ) {
      // Reset
      heatedElements.setCteValue( 0 );
      vector<int> collidedEls = hc.collide( p1, p2 );
      for (int ielem : collidedEls ) {
        if (domain->activeElements[ielem]) {
          heatedElements[ielem] = 1;
        }
      }
      // This is here for the time being
      computePowerDensity();
    }

    void computeHeatedVolume() {
      heatedVolume = 0.0;
      vector<int> indicesHeatedEls = heatedElements.getIndices();
      for (int ielem : indicesHeatedEls) {
        mesh::Element e = domain->mesh->getElement( ielem );
        heatedVolume += e.vol;
      }
    }

    void computePowerDensity() {
      computeHeatedVolume();
      pd = power / heatedVolume;

      //BDEBUG
      elementPulse.setCteValue( 0.0 );
      vector<int> indicesHeatedEls = heatedElements.getIndices();
      for (int ielem : indicesHeatedEls) {
        elementPulse[ielem] = pd;
      }
      //EDEBUG
    }
};
}
#endif
