#include "mesh/Element.h"
#include "Problem.h"
#include "../external/pybind11/include/pybind11/pybind11.h"
#include "../external/pybind11/include/pybind11/stl.h"
#include "../external/pybind11/include/pybind11/eigen.h"
#include "../external/pybind11/include/pybind11/functional.h"

namespace py = pybind11;

PYBIND11_MODULE(MovingHeatSource, m) {
    py::class_<Problem>(m, "Problem", py::dynamic_attr())
        .def(py::init<Problem>())//copy constructor
        .def(py::init<mesh::Mesh&, py::dict&>())
        .def("preIterate", &Problem::preIterate)
        .def("preAssemble", &Problem::preAssemble)
        .def("updateForcedDofs", &Problem::updateForcedDofs)
        .def("assemble", &Problem::assemble)
        .def("iterate", &Problem::iterate)
        .def("gather", &Problem::gather)
        .def("postIterate", &Problem::postIterate)
        .def("initializeIntegrator", &Problem::initializeIntegrator)
        .def("updateFRFpos", &Problem::updateFRFpos)
        .def("setTime", &Problem::setTime)
        .def("setAdvectionSpeed", &Problem::setAdvectionSpeed)
        .def("gather", &Problem::gather)
        .def_readonly("myls", &Problem::myls)
        .def_readonly("ls", &Problem::ls)
        .def_readonly("dofNumbering", &Problem::dofNumbering)
        .def_readonly("freeDofsNumbering", &Problem::freeDofsNumbering)
        .def_readonly("unknown", &Problem::unknown)
        .def_readwrite("previousValues", &Problem::previousValues)
        .def_readonly("mhs", &Problem::mhs)
        .def_readonly("domain", &Problem::domain)
        .def_readwrite("time", &Problem::time)
        .def_readwrite("hasPreIterated", &Problem::hasPreIterated)
        .def_readonly("dt", &Problem::dt)
        .def_readonly("isAdvection", &Problem::isAdvection)
        .def_readonly("advectionSpeed", &Problem::advectionSpeed)
        .def_readonly("dirichletNodes", &Problem::dirichletNodes)
        .def_readonly("gammaNodes", &Problem::gammaNodes)
        .def("interpolate2dirichlet", &Problem::interpolate2dirichlet)
        .def("setAssembling2External", &Problem::setAssembling2External)
        .def("setDeltaT", &Problem::setDeltaT)
        .def("setPointers", &Problem::setPointers)
        .def("setStabilization", &Problem::setStabilization)
        .def("setDirichlet", static_cast<void (Problem::*)(vector<int>, std::function<double(Eigen::Vector3d)>)>(&Problem::setDirichlet),
            "Set Dirichlet condition from indices of facets and function.")
        .def("setDirichlet", static_cast<void (Problem::*)(const vector<int>&, const vector<double>&)>(&Problem::setDirichlet),
            "Set Dirichlet condition from indices of facets and function.")
        .def("setNeumann", static_cast<void (Problem::*)(vector<vector<int>>, double)>(&Problem::setNeumann),
            "Set Neumann condition from array of nodes.")
        .def("setNeumann", static_cast<void (Problem::*)(Eigen::Vector3d, Eigen::Vector3d, double)>(&Problem::setNeumann),
            "Set Neumann condition from plane.")
        .def("setNeumann", static_cast<void (Problem::*)(vector<int>, std::function<Eigen::Vector3d(Eigen::Vector3d)>)>(&Problem::setNeumann),
            "Set Neumann condition from index of facet and flux function.")
        .def("setGamma2Dirichlet", &Problem::setGamma2Dirichlet)
        .def("assembleNeumannGamma", &Problem::assembleNeumannGamma)
        .def("assembleDirichletGamma", &Problem::assembleDirichletGamma)
        .def("clearBCs", &Problem::clearBCs)
        .def("project", &Problem::project)
        .def("getActiveInExternal", static_cast<mesh::MeshTag<int> (Problem::*)( const Problem &, double)>(&Problem::getActiveInExternal),
            "Find interface between two problems")
        .def("findGamma", static_cast<void (Problem::*)( const Problem &)>(&Problem::findGamma),
            "Find interface between two problems")
        .def("findGamma", static_cast<void (Problem::*)( mesh::MeshTag<int>& )>(&Problem::findGamma),
            "Find interface between two problems given external activation MeshTag")
        .def("substractExternal", static_cast<void (Problem::*)( const Problem &, bool, bool)>(&Problem::substractExternal),
            "Substract external domain from domain.")
        .def("intersectExternal", static_cast<void (Problem::*)( const Problem &, bool, bool)>(&Problem::intersectExternal),
            "Intersect external domain with domain.");
    py::class_<mesh::ActiveMesh>(m, "ActiveMesh")
        .def(py::init<mesh::Mesh*>())
        .def("dim", &mesh::ActiveMesh::dim)
        .def_readonly("activeNodes", &mesh::ActiveMesh::activeNodes)
        .def_readonly("activeElements", &mesh::ActiveMesh::activeElements)
        .def_readonly("justActivatedBoundary", &mesh::ActiveMesh::justActivatedBoundary)
        .def_readonly("boundaryFacets", &mesh::ActiveMesh::boundaryFacets)
        .def_readonly("mesh", &mesh::ActiveMesh::mesh)
        .def("setActivation", &mesh::ActiveMesh::setActivation)
        .def("findOwnerElement", &mesh::ActiveMesh::findOwnerElement);
    py::class_<mesh::MeshTag<int>>(m, "MeshTag")//TODO: do it in a loop
        .def(py::init<const mesh::Mesh*>())
        .def(py::init<const mesh::Mesh*, int>())
        .def(py::init<const mesh::Mesh*, int, vector<int>>())
        .def(py::init<const mesh::Mesh*, const std::vector<int>&, const std::vector<int> &, const int>())
        .def("dim", &mesh::MeshTag<int>::dim)
        .def("setValues", &mesh::MeshTag<int>::setValues)
        .def("getIndices", &mesh::MeshTag<int>::getIndices)
        .def_readonly("x", &mesh::MeshTag<int>::x);
    py::class_<mesh::Mesh>(m, "Mesh", py::dynamic_attr())
        .def(py::init<const mesh::Mesh&>())
        .def(py::init<const py::dict&>())
        .def_readonly("pos", &mesh::Mesh::pos)
        .def_readonly("posFRF", &mesh::Mesh::posFRF)
        .def_readonly("con_CellPoint", &mesh::Mesh::con_CellPoint)
        .def_readonly("con_CellCell", &mesh::Mesh::con_CellCell)
        .def_readonly("con_FacetPoint", &mesh::Mesh::con_FacetPoint)//Debugging
        .def_readonly("con_FacetCell", &mesh::Mesh::con_FacetCell)//Debugging
        .def_readonly("nels", &mesh::Mesh::nels)
        .def_readonly("nnodes", &mesh::Mesh::nnodes)
        .def_readonly("speedFRF", &mesh::Mesh::speedFRF)
        .def_readonly("shiftFRF", &mesh::Mesh::shiftFRF)
        .def_readonly("dim", &mesh::Mesh::dim)
        .def_readonly("elementTypes", &mesh::Mesh::elementTypes)
        .def("setSpeedFRF", &mesh::Mesh::setSpeedFRF)
        .def("findOwnerElement", &mesh::Mesh::findOwnerElement)
        .def("getElement", &mesh::Mesh::getElement);
    py::class_<fem::Function>(m, "Function", py::dynamic_attr())
        .def(py::init<const mesh::ActiveMesh*>())
        .def(py::init<const mesh::ActiveMesh*, const Eigen::VectorXd&>())
        .def("evaluate", &fem::Function::evaluate)
        .def("evaluateGrad", &fem::Function::evaluateGrad)
        .def("interpolate", &fem::Function::interpolate)
        .def("setValues", &fem::Function::setValues)
        .def_readonly("values", &fem::Function::values);
    //This export won't work unless list<Function> is made into
    //an opaque type or interpolate is wrapped into something else
    //m.def( "interpolate", &fem::interpolate, "interpolate list of sourceFunctions to targetFunctions" );
    py::class_<mesh::Connectivity>(m, "Connectivity", py::dynamic_attr())
        .def("getLocalCon", &mesh::Connectivity::getLocalCon)
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
        .def("getCentroid", &mesh::Element::getCentroid)
        .def("getSizeAlongVector", &mesh::Element::getSizeAlongVector);
    py::enum_<ElementType>(m, "ElementType")
        .value("point1", ElementType::point1)
        .value("line2", ElementType::line2)
        .value("triangle3", ElementType::triangle3)
        .value("quad4", ElementType::quad4);
    py::class_<HeatSource>(m, "HeatSource")
        .def(py::init<>())
        .def_readwrite("currentPosition", &HeatSource::currentPosition)
        .def_readonly("pulse", &HeatSource::pulse)
        .def_readonly("speed", &HeatSource::speed)
        .def("setSpeed", &HeatSource::setSpeed);
    py::class_<LinearSystem>(m, "LinearSystem")
        .def(py::init<>())
        .def_readonly("lhs", &LinearSystem::lhs)
        .def_readonly("rhs", &LinearSystem::rhs)
        .def_readonly("sol", &LinearSystem::sol)
        .def("assemble", &LinearSystem::assemble)
        .def("solve", &LinearSystem::solve)
        .def("cleanup", &LinearSystem::cleanup)
        .def("ndofs", &LinearSystem::getNdofs)
        .def( py::init<Problem&, Problem&>() );
}
