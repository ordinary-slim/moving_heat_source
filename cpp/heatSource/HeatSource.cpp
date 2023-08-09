#include "HeatSource.h"
#include "../Problem.h"

namespace py = pybind11;

heat::HeatSource::HeatSource( pybind11::dict &input, Problem *problem ) {
  this->radius = py::cast<double>(input["radius"]);
  this->power = py::cast<double>(input["power"]);
  this->speed = CreateEigenVector(py::array_t<double>(input["HeatSourceSpeed"]));
  this->position = CreateEigenVector(py::array_t<double>(input["initialPosition"]));

  if (input.contains("efficiency")) this->efficiency = py::cast<double>(input["efficiency"]);

  this->problem = problem;
  this->pulse.resize( problem->domain.mesh->nnodes );
}

void heat::HeatSource::preIterate() {
  if (path != NULL) {
    currentTrack = path->interpolateTrack( problem->time );
    if (currentTrack != NULL) {
      this->speed = currentTrack->getSpeed();
      this->power = currentTrack->power;
    } else {
      throw std::invalid_argument("Time is out of bounds.");
    }
  }
  step( problem->dt );
}
