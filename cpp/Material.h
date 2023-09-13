#include <pybind11/stl.h>
#include <iostream>

namespace py = pybind11;

struct ThermalMaterial {
  double density = 1.0;
  double conductivity = 1.0;
  double specificHeat = 1.0;
  double convectionCoeff = 0.0;

  ThermalMaterial( py::dict &input ) {
    density = py::cast<double>(input["rho"]);
    conductivity = py::cast<double>(input["conductivity"]);
    specificHeat = py::cast<double>(input["specific_heat"]);
    if (input.contains("convectionCoeff")) {
      convectionCoeff = py::cast<double>(input["convectionCoeff"]);
    }
  }
};
