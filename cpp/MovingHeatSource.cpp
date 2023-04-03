#include "mesh/Element.h"
#include "Problem.h"
#include "../external/pybind11/include/pybind11/pybind11.h"
#include "../external/pybind11/include/pybind11/stl.h"
#include "../external/pybind11/include/pybind11/eigen.h"

namespace py = pybind11;

PYBIND11_MODULE(MovingHeatSource, m) {
    py::class_<Problem>(m, "Problem", py::dynamic_attr())
        .def(py::init<>())
        .def(py::init<Problem>())
        .def("initialize", &Problem::initialize)
        .def("initializeIntegrator", &Problem::initializeIntegrator)
        .def("updateFRFpos", &Problem::updateFRFpos)
        .def("iterate", &Problem::iterate)
        .def("preIterate", &Problem::preIterate)
        .def("postIterate", &Problem::postIterate)
        .def("setTime", &Problem::setTime)
        .def("setAdvectionSpeed", &Problem::setAdvectionSpeed)
        .def_readonly("unknown", &Problem::unknown)
        .def_readonly("mhs", &Problem::mhs)
        .def_readonly("mesh", &Problem::mesh)
        .def_readwrite("time", &Problem::time)
        .def_readonly("dt", &Problem::dt)
        .def_readonly("isAdvection", &Problem::isAdvection)
        .def_readonly("advectionSpeed", &Problem::advectionSpeed)
        .def("setPointers", &Problem::setPointers)
        .def("setStabilization", &Problem::setStabilization)
        .def("setNeumann", static_cast<void (Problem::*)(vector<vector<int>>, double)>(&Problem::setNeumann),
            "Set Neumann condition from array of nodes.")
        .def("setNeumann", static_cast<void (Problem::*)(Eigen::Vector3d, Eigen::Vector3d, double)>(&Problem::setNeumann),
            "Set Neumann condition from plane.")
        .def("activateDomain", &Problem::activateDomain);
    py::class_<mesh::Mesh>(m, "Mesh", py::dynamic_attr())
        .def(py::init<>())
        .def_readonly("pos", &mesh::Mesh::pos)
        .def_readonly("posFRF", &mesh::Mesh::posFRF)
        .def_readonly("con_CellPoint", &mesh::Mesh::con_CellPoint)
        .def_readonly("con_CellCell", &mesh::Mesh::con_CellCell)
        .def_readonly("con_FacetPoint", &mesh::Mesh::con_FacetPoint)//Debugging
        .def_readonly("con_FacetCell", &mesh::Mesh::con_FacetCell)//Debugging
        .def_readonly("nels", &mesh::Mesh::nels)
        .def_readonly("nnodes", &mesh::Mesh::nnodes)
        .def_readonly("activeNodes", &mesh::Mesh::activeNodes)
        .def_readonly("activeElements", &mesh::Mesh::activeElements)
        .def_readonly("shiftFRF", &mesh::Mesh::shiftFRF)
        .def("setSpeedFRF", &mesh::Mesh::setSpeedFRF)
        .def("findOwnerElement", &mesh::Mesh::findOwnerElement)
        .def("getElement", &mesh::Mesh::getElement);
    py::class_<FEMFunction>(m, "FEMFunction", py::dynamic_attr())
        .def(py::init<>())
        .def("evaluateVal", &FEMFunction::evaluateVal)
        .def("evaluateValNPrevVals", &FEMFunction::evaluateValNPrevVals)
        .def("getFromExternal", &FEMFunction::getFromExternal)
        .def("forceFromExternal", &FEMFunction::forceFromExternal)
        .def("releaseDirichlet", &FEMFunction::releaseDirichlet)
        .def_readonly("values", &FEMFunction::values);
    py::class_<mesh::Connectivity>(m, "Connectivity", py::dynamic_attr())
        .def_readonly("con", &mesh::Connectivity::con)
        .def_readonly("nels_oDim", &mesh::Connectivity::nels_oDim)
        .def_readonly("nels_tDim", &mesh::Connectivity::nels_tDim);
    m.def( "transpose", &mesh::transpose, "(d -> d') ---> (d' -> d)" );
    m.def( "intersect", &mesh::intersect, "(d -> d''), (d'' -> d') ---> (d -> d')" );
    m.def( "build", &mesh::build, "(D -> 0), (D -> D) ---> (d -> 0), (D -> d)" );
    py::class_<mesh::Element>(m, "mesh::Element", py::dynamic_attr())
        .def(py::init<>())
        .def_readonly("pos", &mesh::Element::pos)
        .def_readonly("nnodes", &mesh::Element::nnodes)
        .def_readonly("gpos", &mesh::Element::gpos)
        .def_readonly("ngpoints", &mesh::Element::ngpoints)
        .def_readonly("gpweight", &mesh::Element::gpweight)
        .def_readonly("con", &mesh::Element::con)
        .def_readonly("vol", &mesh::Element::vol)
        .def_readonly("dimension", &mesh::Element::dim)
        .def_readonly("elementType", &mesh::Element::elementType)
        .def_readonly("GradBaseGpVals", &mesh::Element::GradBaseGpVals)
        .def("computeDerivatives", &mesh::Element::computeDerivatives)
        .def("getSizeAlongVector", &mesh::Element::getSizeAlongVector);
    py::class_<HeatSource>(m, "HeatSource")
        .def(py::init<>())
        .def_readwrite("currentPosition", &HeatSource::currentPosition)
        .def_readonly("pulse", &HeatSource::pulse)
        .def_readonly("speed", &HeatSource::speed)
        .def("setSpeed", &HeatSource::setSpeed);
}
