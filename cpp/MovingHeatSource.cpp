#include "mesh/Element.h"
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
        .def_readonly("unknown", &Problem::unknown)
        .def_readonly("mhs", &Problem::mhs)
        .def_readonly("mesh", &Problem::mesh)
        .def_readwrite("time", &Problem::time)
        .def_readonly("dt", &Problem::dt)
        .def_readonly("isAdvection", &Problem::isAdvection)
        .def_readonly("advectionSpeed", &Problem::advectionSpeed)
        .def("activateDomain", &Problem::activateDomain);
    py::class_<mesh::Mesh>(m, "Mesh", py::dynamic_attr())
        .def(py::init<>())
        .def_readonly("pos", &mesh::Mesh::pos)
        .def_readonly("pos_noAdv", &mesh::Mesh::pos_noAdv)
        .def_readonly("con_CellPoint", &mesh::Mesh::con_CellPoint)
        .def_readonly("con_CellCell", &mesh::Mesh::con_CellCell)
        .def_readonly("con_FacetPoint", &mesh::Mesh::con_FacetPoint)//Debugging
        .def_readonly("con_FacetCell", &mesh::Mesh::con_FacetCell)//Debugging
        .def_readonly("nels", &mesh::Mesh::nels)
        .def_readonly("nnodes", &mesh::Mesh::nnodes)
        .def_readonly("activeElements", &mesh::Mesh::activeElements)
        .def("generate1DMesh", &mesh::Mesh::generate1DMesh)
        .def("findOwnerElement", &mesh::Mesh::findOwnerElement)
        .def("getElement", &mesh::Mesh::getElement);
    py::class_<FEMFunction>(m, "FEMFunction", py::dynamic_attr())
        .def(py::init<>())
        .def("interpolate", &FEMFunction::interpolate)
        .def_readonly("values", &FEMFunction::values);
    py::class_<mesh::Connectivity>(m, "Connectivity", py::dynamic_attr())
        .def_readonly("con", &mesh::Connectivity::con)
        .def_readonly("nels_oDim", &mesh::Connectivity::nels_oDim)
        .def_readonly("nels_tDim", &mesh::Connectivity::nels_tDim);
    m.def( "transpose", &mesh::transpose, "(d -> d') ---> (d' -> d)" );
    m.def( "intersect", &mesh::intersect, "(d -> d''), (d'' -> d') ---> (d -> d')" );
    m.def( "build", &mesh::build, "(D -> 0), (D -> D) ---> (d -> 0), (D -> d)" );
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
