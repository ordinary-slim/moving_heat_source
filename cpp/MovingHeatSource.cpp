#include "mesh/Element.h"
#include "Problem.h"
#include "additiveManufacturing/Printer.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/operators.h>
#include "linearAlgebra/LinearSystem.h"
#include "linearAlgebra/Solver.h"
#ifdef hasMATLAB
#include "linearAlgebra/MATLAB.h"
#endif

namespace py = pybind11;

PYBIND11_MODULE(cpp, m) {
    py::class_<Problem>(m, "Problem", py::dynamic_attr())
        //.def(py::init<Problem>())//copy constructor
        // doesnt work currently with HeatSource unique pointer
        .def(py::init<mesh::Mesh&, py::dict&>())
        .def("preIterate", &Problem::preIterate, "Time goes from t^n to t^{n+1}",
            py::arg("canPreassemble") = true)
        .def("preAssemble", &Problem::preAssemble, "Before assembly operations",
            py::arg("allocateLs") = true )
        .def("updateForcedDofs", &Problem::updateForcedDofs)
        .def("assemble", &Problem::assemble,
            py::arg("externalProblem") = static_cast<Problem *>(nullptr))
        .def("iterate", &Problem::iterate)
        .def("gather", &Problem::gather)
        .def("postIterate", &Problem::postIterate)
        .def("initializeIntegrator", &Problem::initializeIntegrator)
        .def("setAdvectionSpeed", &Problem::setAdvectionSpeed)
        .def("gather", &Problem::gather)
        .def_readonly("stabilizationScheme", &Problem::stabilizationScheme)
        .def_readonly("isCoupled", &Problem::isCoupled)
        .def_readonly("ls", &Problem::ls)
        .def_readonly("dofNumbering", &Problem::dofNumbering)
        .def_readonly("freeDofsNumbering", &Problem::freeDofsNumbering)
        .def_readonly("unknown", &Problem::unknown)
        .def_readwrite("previousValues", &Problem::previousValues)
        .def_property_readonly("mhs", [](const Problem& p){ return p.mhs.get(); },
            py::return_value_policy::reference_internal)
        .def_readonly("domain", &Problem::domain)
        .def_readonly("materials", &Problem::materials)
        .def_readwrite("time", &Problem::time)
        .def_readwrite("hasPreIterated", &Problem::hasPreIterated)
        .def_readonly("dt", &Problem::dt)
        .def_readonly("isAdvection", &Problem::isAdvection)
        .def_readonly("advectionSpeed", &Problem::advectionSpeed)
        .def_readonly("dirichletNodes", &Problem::dirichletNodes)
        .def_readonly("gammaNodes", &Problem::gammaNodes)
        .def_readonly("forcedDofs", &Problem::forcedDofs)
        .def("interpolate2dirichlet", &Problem::interpolate2dirichlet)
        .def("setDt", &Problem::setDt)
        .def("setPointers", &Problem::setPointers)
        .def("setStabilization", &Problem::setStabilization)
        .def("setDirichlet", static_cast<void (Problem::*)(vector<int>, std::function<double(Eigen::Vector3d)>)>(&Problem::setDirichlet),
            "Set Dirichlet condition from indices of facets and function.")
        .def("setDirichlet", static_cast<void (Problem::*)(const vector<int>&, const vector<double>&)>(&Problem::setDirichlet),
            "Set Dirichlet condition from indices of facets and function.")
        .def("setNeumann", static_cast<void (Problem::*)(vector<vector<unsigned int>>, double)>(&Problem::setNeumann),
            "Set Neumann condition from array of nodes.")
        .def("setNeumann", static_cast<void (Problem::*)(Eigen::Vector3d, Eigen::Vector3d, double)>(&Problem::setNeumann),
            "Set Neumann condition from plane.")
        .def("setNeumann", static_cast<void (Problem::*)(vector<int>, std::function<Eigen::Vector3d(Eigen::Vector3d)>)>(&Problem::setNeumann),
            "Set Neumann condition from index of facet and flux function.")
        .def("setConvection", &Problem::setConvection,
            "Set all external boundary to convection.",
            py::arg("resetBcs") = true)
        .def("setGamma2Neumann", &Problem::setGamma2Neumann)
        .def("setGamma2Dirichlet", &Problem::setGamma2Dirichlet)
        .def("clearBCs", &Problem::clearBCs)
        .def("clearGamma", &Problem::clearGamma)
        .def("getActiveInExternal", &Problem::getActiveInExternal, py::return_value_policy::take_ownership)
        .def("updateInterface", static_cast<void (Problem::*)( const Problem &)>(&Problem::updateInterface),
            "Find interface between two problems")
        .def("updateInterface", static_cast<void (Problem::*)( mesh::MeshTag<int>& )>(&Problem::updateInterface),
            "Find interface between two problems given external activation MeshTag")
        .def("substractExternal", &Problem::substractExternal,
            "Substract external domain from domain.",
            py::arg("pExt"), py::arg("updateGamma") = true )
        .def("intersectExternal", &Problem::intersectExternal,
            "Intersect domain with external domain.",
            py::arg("pExt"), py::arg("updateGamma") = true, py::arg("interpolateIntersected") = false )
        .def("uniteExternal", &Problem::uniteExternal,
            "Unite domain with external domain.",
            py::arg("pExt"), py::arg("updateGamma") = true );
    py::class_<LinearSystem, std::shared_ptr<LinearSystem>>(m, "LinearSystem")
        .def( py::init<>() )
        .def( py::init<Problem&>() )
        .def( py::init<Problem&, Problem&>() )
        .def_static("Create", &LinearSystem::Create, py::return_value_policy::reference)
        .def_readonly("lhs", &LinearSystem::lhs)
        .def_readonly("rhs", &LinearSystem::rhs)
        .def_readonly("sol", &LinearSystem::sol)
        .def("setInitialGuess", &LinearSystem::setInitialGuess,
            py::arg("p1"), py::arg("p2") = static_cast<Problem *>(nullptr))
        .def("assemble", &LinearSystem::assemble)
        .def_property_readonly("solver", [](const LinearSystem& ls){ return ls.solver.get(); },
            py::return_value_policy::reference_internal)
        .def("solve", &LinearSystem::solve)
        .def("cleanup", &LinearSystem::cleanup)
        .def("setSolver", &LinearSystem::setSolver)
        .def("ndofs", &LinearSystem::getNdofs);
    py::class_<Solver>(m, "Solver");
#ifdef hasMATLAB
    py::class_<MatLabSession>(m, "MatLabSession")
        .def( py::init<>() );
    py::class_<MatLabSolver>(m, "MatLabSolver")
        .def("setSession", &MatLabSolver::setSession, py::arg("matlabSession") = static_cast<MatLabSession *>(nullptr));
#endif
    py::class_<mesh::Domain>(m, "Domain")
        .def_readonly("posLab", &mesh::Domain::posLab)
        .def_readonly("translationLab", &mesh::Domain::translationLab)
        .def_readonly("speedDomain", &mesh::Domain::speedDomain)
        .def("setSpeed", &mesh::Domain::setSpeed)
        .def("dim", &mesh::Domain::dim)
        .def_readonly("materialTag", &mesh::Domain::materialTag)
        .def_readonly("activeNodes", &mesh::Domain::activeNodes)
        .def_readonly("activeElements", &mesh::Domain::activeElements)
        .def_readonly("justActivatedBoundary", &mesh::Domain::justActivatedBoundary)
        .def_readonly("boundaryFacets", &mesh::Domain::boundaryFacets)
        .def_readonly("mesh", &mesh::Domain::mesh)
        .def("assembleMassMatrix", &mesh::Domain::assembleMassMatrix)
        .def("project", &mesh::Domain::project)
        .def("setMaterialSets", &mesh::Domain::setMaterialSets)
        .def("setActivation", &mesh::Domain::setActivation)
        .def("resetActivation", &mesh::Domain::resetActivation)
        .def("deactivate", &mesh::Domain::deactivate)
        .def("intersect", &mesh::Domain::intersect)
        .def("inPlaneRotate", &mesh::Domain::inPlaneRotate)
        .def("projectCellTag", &mesh::Domain::projectCellTag)
        .def("findOwnerElements", &mesh::Domain::findOwnerElements);
    py::class_<mesh::MeshTag<int>>(m, "MeshTag")//TODO: do it in a loop
        .def(py::init<const mesh::MeshTag<int>&>())
        .def(py::init<const mesh::Mesh*>())
        .def(py::init<const mesh::Mesh*, int>())
        .def(py::init<const mesh::Mesh*, int, vector<int>>())
        .def("__not__", [](const mesh::MeshTag<int> &mt) {
            return !mt;
        }, py::is_operator())
        .def("__getitem__", [](mesh::MeshTag<int>& t, size_t index) { return &t[index];})
        .def("__len__", [](mesh::MeshTag<int>& t) { return t.size();})
        .def("dim", &mesh::MeshTag<int>::dim)
        .def("setIndices", &mesh::MeshTag<int>::setIndices)
        .def("getIndices", &mesh::MeshTag<int>::getIndices)
        .def_readonly("x", &mesh::MeshTag<int>::x);
    py::class_<ThermalMaterial>(m, "ThermalMaterial")//TODO: do it in a loop
        .def_readonly("conductivity", &ThermalMaterial::conductivity)
        .def_readonly("density", &ThermalMaterial::density)
        .def_readonly("specificHeat", &ThermalMaterial::specificHeat);
    py::class_<mesh::Mesh>(m, "Mesh", py::dynamic_attr())
        //.def(py::init<const mesh::Mesh&>()) AABB_tree doesnt allow this
        .def(py::init<const py::dict&>())
        .def_readonly("pos", &mesh::Mesh::pos)
        .def_readonly("con_CellPoint", &mesh::Mesh::con_CellPoint)
        .def_readonly("con_CellCell", &mesh::Mesh::con_CellCell)
        .def_readonly("con_FacetPoint", &mesh::Mesh::con_FacetPoint)//Debugging
        .def_readonly("con_FacetCell", &mesh::Mesh::con_FacetCell)//Debugging
        .def_readonly("nels", &mesh::Mesh::nels)
        .def_readonly("nnodes", &mesh::Mesh::nnodes)
        .def_readonly("dim", &mesh::Mesh::dim)
        .def("getElementType", &mesh::Mesh::getElementType,
            py::arg("ielem") = 0)
        .def("findOwnerElements", &mesh::Mesh::findOwnerElements)
        .def("findCollidingElements", static_cast<std::vector<int> (mesh::Mesh::*)(const MyAABB&) const>(&mesh::Mesh::findCollidingElements))
        .def("findCollidingElements", static_cast<std::vector<int> (mesh::Mesh::*)(const MyOBB&) const>(&mesh::Mesh::findCollidingElements))
        .def("findCollidingElements", static_cast<std::vector<int> (mesh::Mesh::*)(const Eigen::Vector3d&, const double R) const>(&mesh::Mesh::findCollidingElements))
        .def("getElementGeometry", &mesh::Mesh::getElementGeometry)
        .def("getElement", &mesh::Mesh::getElement);
    py::class_<AbstractFunction>(m, "AbstractFunction");
    py::class_<fem::Function, AbstractFunction>(m, "Function")
        .def(py::init<const mesh::Domain*>())
        .def(py::init<const mesh::Domain*, Eigen::VectorXd&>())
        .def(py::init<const fem::Function&>())
        .def(py::self - py::self)
        .def(py::self / double())
        .def("evaluate", &fem::Function::evaluate)
        .def("evaluateGrad", &fem::Function::evaluateGrad)
        .def("interpolate", py::overload_cast<const AbstractFunction &, bool>(&fem::Function::interpolate),
          py::arg("extFEMFunc"), py::arg("ignoreOutside") = false )
        .def("interpolate", py::overload_cast<const AbstractFunction &, const mesh::MeshTag<int> &, bool>(&fem::Function::interpolate),
          py::arg("extFEMFunc"), py::arg("nodalTag"), py::arg("ignoreOutside") = false )
        .def("interpolateInactive", &fem::Function::interpolateInactive,
          py::arg("extFEMFunc"), py::arg("ignoreOutside") = false )
        .def("getL2Norm", &fem::Function::getL2Norm)
        .def_readonly("values", &fem::Function::values);
    m.def( "interpolate", &fem::interpolate);
    py::class_<mesh::Connectivity>(m, "Connectivity", py::dynamic_attr())
        .def("getLocalCon", &mesh::Connectivity::getLocalCon)
        .def_readonly("con", &mesh::Connectivity::con)
        .def_readonly("nels_oDim", &mesh::Connectivity::nels_oDim)
        .def_readonly("nels_tDim", &mesh::Connectivity::nels_tDim);
    m.def( "transpose", &mesh::transpose, "(d -> d') ---> (d' -> d)" );
    m.def( "intersect", &mesh::intersect, "(d -> d''), (d'' -> d') ---> (d -> d')" );
    m.def( "buildBoundaryConnectivities", &mesh::buildBoundaryConnectivities, "(D -> 0), (D -> D) ---> (D-1 -> 0), (D -> D-1)" );
    py::class_<mesh::Element>(m, "Element", py::dynamic_attr())
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
        .def("getCentroid", &mesh::Element::getCentroid)
        .def("getSizeAlongVector", &mesh::Element::getSizeAlongVector)
        .def("evaluateShaFuns", &mesh::Element::evaluateShaFuns)
        .def("evaluateGradShaFuns", &mesh::Element::evaluateGradShaFuns);
    py::enum_<ElementType>(m, "ElementType")
        .value("point1", ElementType::point1)
        .value("line2", ElementType::line2)
        .value("triangle3", ElementType::triangle3)
        .value("quad4", ElementType::quad4)
        .value("hexa8", ElementType::hexa8);
    py::class_<MyAABB>(m, "MyAABB")
        .def(py::init<Eigen::Vector3d, Eigen::Vector3d, bool>());
    py::class_<MyOBB>(m, "MyOBB")
        .def(py::init<Eigen::Vector3d, Eigen::Vector3d, double, double>());
    py::enum_<heat::HeatSourceType>(m, "HeatSourceType")
        .value("none", heat::HeatSourceType::none)
        .value("gaussian1d", heat::HeatSourceType::gaussian1d)
        .value("gaussian2d", heat::HeatSourceType::gaussian2d)
        .value("gaussian3d", heat::HeatSourceType::gaussian3d)
        .value("constant", heat::HeatSourceType::constant)
        .value("lumped", heat::HeatSourceType::lumped);
    py::class_<heat::HeatSource>(m, "HeatSource")
        .def_readonly("position", &heat::HeatSource::position)
        .def_readonly("pulse", &heat::HeatSource::pulse)
        .def_readonly("power", &heat::HeatSource::power)
        .def_readonly("speed", &heat::HeatSource::speed)
        .def_readonly("radius", &heat::HeatSource::radius)
        .def_readonly("currentTrack", &heat::HeatSource::currentTrack)
        .def_readonly("type", &heat::HeatSource::type)
        .def_property_readonly("path", [](const heat::HeatSource& h){ return h.path.get(); },
            py::return_value_policy::reference_internal)
        .def("__call__", &heat::HeatSource::operator())
        .def("setPower", &heat::HeatSource::setPower)
        .def("setPosition", &heat::HeatSource::setPosition)
        .def("setSpeed", &heat::HeatSource::setSpeed)
        .def("setPath", &heat::HeatSource::setPath)
        .def("getNextTrack", &heat::HeatSource::getNextTrack, pybind11::return_value_policy::reference);
    py::class_<heat::LumpedHeatSource, heat::HeatSource>(m, "LumpedHeatSource")
        .def_readonly("heatedElements", &heat::LumpedHeatSource::heatedElements)
        .def_readonly("elementPulse", &heat::LumpedHeatSource::elementPulse)//DEBUG
        .def("markHeatedElements", &heat::LumpedHeatSource::markHeatedElements);
    py::class_<HatchCollider>(m, "HatchCollider");
    py::class_<Printer, HatchCollider>(m, "Printer")
        .def( py::init<Problem*, double, double, double>(),
            py::arg("problem"), py::arg("width"), py::arg("height"), py::arg("depth")=0.0 )
        .def("collide", &Printer::collide)
        .def("setDepositionTemperature", &Printer::setDepositionTemperature)
        .def("melt", &Printer::melt,
            py::arg("p1"), py::arg("p2"), py::arg("materialTag") = static_cast<mesh::MeshTag<int> *>(nullptr),
            py::arg("idxTargetSet") = 0)
        .def("deposit", &Printer::deposit,
            py::arg("p1"), py::arg("p2"), py::arg("activeEls") = static_cast<mesh::MeshTag<int> *>(nullptr),
            py::arg("modifyValues") = true);
    py::enum_<heat::TrackType>(m, "TrackType")
        .value("printing", heat::TrackType::printing)
        .value("cooling", heat::TrackType::cooling)
        .value("dwelling", heat::TrackType::dwelling)
        .value("recoating", heat::TrackType::recoating);
    py::class_<heat::Track>(m, "Track")
        .def_readonly("p0", &heat::Track::p0)
        .def_readonly("p1", &heat::Track::p1)
        .def("getSpeed", &heat::Track::getSpeed)
        .def("isOver", &heat::Track::isOver)
        .def_readonly("speed", &heat::Track::speed)
        .def_readonly("power", &heat::Track::power)
        .def_readonly("type", &heat::Track::type)
        .def_readonly("isNewX", &heat::Track::isNewX)
        .def_readonly("isNewY", &heat::Track::isNewY)
        .def_readonly("isNewZ", &heat::Track::isNewZ)
        .def_readonly("startTime", &heat::Track::startTime)
        .def_readonly("endTime", &heat::Track::endTime);
    py::class_<heat::Path>(m, "Path")
        .def("interpolateTrack", &heat::Path::interpolateTrack, pybind11::return_value_policy::reference)
        .def("interpolatePosition", &heat::Path::interpolatePosition)
        .def("isOver", &heat::Path::isOver);
}
