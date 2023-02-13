#include "domain/element.h"
#include "problem.h"
#include "../external/pybind11/include/pybind11/pybind11.h"
#include "../external/pybind11/include/pybind11/stl.h"
#include "../external/pybind11/include/pybind11/eigen.h"

namespace py = pybind11;

PYBIND11_MODULE(MovingHeatSource, m) {
    py::class_<Problem>(m, "Problem", py::dynamic_attr())
        .def(py::init<>())
        .def("initialize", &Problem::initialize)
        .def("initializeIntegrator", &Problem::initializeIntegrator)
        .def("updateFRF_positions", &Problem::updateFRF_positions)
        .def("iterate", &Problem::iterate)
        .def("setTime", &Problem::setTime)
        .def("setAdvectionSpeed", &Problem::setAdvectionSpeed)
        .def_readonly("solution", &Problem::solution)
        //.def_readwrite("pulse", &Problem::pulse)//debugging
        .def_readonly("mhs", &Problem::mhs)
        .def_readonly("mesh", &Problem::mesh)
        .def_readwrite("time", &Problem::time)
        .def_readonly("dt", &Problem::dt)
        .def_readonly("isAdvection", &Problem::isAdvection)
        .def_readonly("advectionSpeed", &Problem::advectionSpeed)
        .def("activateDomain", &Problem::activateDomain);
    py::class_<Mesh>(m, "Mesh", py::dynamic_attr())
        .def(py::init<>())
        .def_readonly("pos", &Mesh::pos)
        .def_readonly("pos_noAdv", &Mesh::pos_noAdv)
        .def_readonly("con_CellPoint", &Mesh::con_CellPoint)
        .def_readonly("nels", &Mesh::nels)
        .def_readonly("nnodes", &Mesh::nnodes)
        .def_readonly("activeElements", &Mesh::activeElements)
        .def("generate1DMesh", &Mesh::generate1DMesh)
        .def("getElement", &Mesh::getElement);
    py::class_<Connectivity>(m, "Connectivity", py::dynamic_attr())
        .def_readonly("con", &Connectivity::con);
    py::class_<Element>(m, "Element", py::dynamic_attr())
        .def(py::init<>())
        .def_readonly("pos", &Element::pos)
        .def_readonly("nnodes", &Element::nnodes)
        .def_readonly("gpos", &Element::gpos)
        .def_readonly("ngpoints", &Element::ngpoints)
        .def_readonly("gpweight", &Element::gpweight)
        .def_readonly("con", &Element::con)
        .def_readonly("vol", &Element::vol)
        .def_readonly("dimension", &Element::dim)
        .def_readonly("elementType", &Element::elementType)
        .def_readonly("GradBaseGpVals", &Element::GradBaseGpVals)
        .def("computeDerivatives", &Element::computeDerivatives);
    py::class_<HeatSource>(m, "HeatSource")
        .def(py::init<>())
        .def_readwrite("currentPosition", &HeatSource::currentPosition)
        .def_readonly("speed", &HeatSource::speed)
        .def("setSpeed", &HeatSource::setSpeed);
}
