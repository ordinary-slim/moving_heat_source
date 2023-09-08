#include "HeatSource.h"
#include "../Problem.h"

namespace py = pybind11;

namespace heat
{
HeatSource::HeatSource( pybind11::dict &input, Problem *problem ) {
  this->radius = py::cast<double>(input["radius"]);
  this->power = py::cast<double>(input["power"]);
  this->speed = CreateEigenVector(py::array_t<double>(input["HeatSourceSpeed"]));
  this->position = CreateEigenVector(py::array_t<double>(input["initialPosition"]));

  if (input.contains("efficiency")) this->efficiency = py::cast<double>(input["efficiency"]);

  this->problem = problem;
  this->pulse.resize( problem->domain.mesh->nnodes );
}

void HeatSource::setPath( std::vector<Eigen::Vector3d> &coordinates,
      std::vector<double> &speeds,
      std::vector<double> &powers,
      std::vector<int> &arePrinting ) {

    path = std::make_unique<heat::Path>( coordinates, speeds, powers, arePrinting );
    currentTrack = &path->tracks[0];
    position = path->interpolatePosition( problem->time );
}

const Track* HeatSource::getNextTrack() const {
  const Track* track = NULL;
  if (currentTrack->index < (path->tracks.size() + 1)) {
    track = &path->tracks[currentTrack->index + 1];
  }
  return track;
}

void HeatSource::preIterate() {
  if (path != NULL) {
    currentTrack = path->interpolateTrack( problem->time );
    cout << "Current track is now heading to \n " << *currentTrack->p1 << endl;
    if (currentTrack != NULL) {
      this->speed = currentTrack->getSpeed();
      this->power = currentTrack->power;
    } else {
      throw std::invalid_argument("Time is out of bounds.");
    }
  }
  step( problem->dt );
}

void HeatSource::postIterate() {
  /*
   * Placeholder
   */
}

}
