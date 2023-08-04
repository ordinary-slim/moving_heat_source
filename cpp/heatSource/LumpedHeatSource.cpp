#include "../Problem.h"
#include "HeatSource.h"
#include "LumpedHeatSource.h"

heat::LumpedHeatSource::LumpedHeatSource( double width, double height,
          py::dict &input, Problem *problem ) :
  HeatSource( input, problem ),
  hc( problem->domain.mesh, width, height ),
  heatedElements( mesh::MeshTag<int>( problem->domain.mesh, problem->domain.mesh->dim, 0 ) ),
  elementPulse( mesh::MeshTag<double>( problem->domain.mesh, problem->domain.mesh->dim, 0.0 ) )//DEBUG
{
  type = lumped;
}

void heat::LumpedHeatSource::markHeatedElements( const Eigen::Vector3d &p1, const Eigen::Vector3d &p2 ) {
  // Reset
  heatedElements.setCteValue( 0 );
  vector<int> collidedEls = hc.collide( p1, p2 );
  for (int ielem : collidedEls ) {
    if (problem->domain.activeElements[ielem]) {
      heatedElements[ielem] = 1;
    }
  }
  // This is here for the time being
  computePowerDensity();
}

void heat::LumpedHeatSource::computeHeatedVolume() {
  heatedVolume = 0.0;
  vector<int> indicesHeatedEls = heatedElements.getIndices();
  for (int ielem : indicesHeatedEls) {
    mesh::Element e = problem->domain.mesh->getElement( ielem );
    heatedVolume += e.vol;
  }
}

void heat::LumpedHeatSource::computePowerDensity() {
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
