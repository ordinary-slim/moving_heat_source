#include <chrono>
#include "Mesh.h"
#include "../../external/pybind11/include/pybind11/pybind11.h"
#include "../../external/pybind11/include/pybind11/eigen.h"
#include "../../external/pybind11/include/pybind11/stl.h"

namespace py = pybind11;

mesh::Mesh::Mesh(const py::dict &input) {
  //SET ELEMENT TYPE
  ElementType cell_type_flag;
  string aux_cell_type = py::cast<string>( input["cell_type"] );
  if        (aux_cell_type == "line2") {
    cell_type_flag = line2;
  } else if (aux_cell_type == "triangle3") {
    cell_type_flag = triangle3;
  } else if (aux_cell_type == "quad4") {
    cell_type_flag = quad4;
  }
  if (input.contains("numberOfGaussPoints")){
    ngpointsCell = py::cast<int>( input["numberOfGaussPoints"] );
  }
  // reference element. no support for mixed meshes yet
  refCellEl = ReferenceElement(cell_type_flag, ngpointsCell);
  ElementType FacetElType = getIncidentElType(cell_type_flag, refCellEl.dim-1);
  refFacetEl = ReferenceElement(FacetElType);

  //READ POINTS
  py::array points = input["points"];
  nnodes = points.shape(0);
  dim    = points.shape(1);
  pos.resize( nnodes, 3 );
  posFRF.resize( nnodes, 3 );
  pos.setZero();
  auto aux_points = points.unchecked<double>();
  for ( int ipoint = 0; ipoint < nnodes; ipoint++) {
    for ( int idim = 0; idim < dim; idim++) {
      pos(ipoint, idim) =  aux_points(ipoint, idim);
    }
  }
  posFRF = pos;
  //READ ELEMENTS
  py::array cells    = input["cells"];
  nels          = cells.shape(0);
  nnodes_per_el = cells.shape(1);
  con_CellPoint.con.resize( nels );
  elementTypes.resize( nels );
  auto aux_cells = cells.unchecked<int>();//Receiving uint here breaks code
  for ( int icell = 0; icell < nels; icell++) {
    con_CellPoint.con[icell].resize( nnodes_per_el );
    for ( int inode = 0; inode < nnodes_per_el; inode++) {
      con_CellPoint.con[icell][inode] =  aux_cells(icell, inode);
    }
  }
  std::fill (elementTypes.begin(), elementTypes.end(), 
      cell_type_flag);
  // reference element. no support for mixed s yet

  //CONNECTIVITIES
  //Manually fill D0 connectivity fields
  //Very unsafe! Need privatization
  con_CellPoint.oDim = refCellEl.dim;
  con_CellPoint.tDim = 0;
  con_CellPoint.nels_oDim = nels;
  con_CellPoint.nels_tDim = nnodes;
  con_CellPoint.oelType = refCellEl.elementType;
  con_CellPoint.telType = point1;
  printf("Nels = %i\n", nels);
  cout << "Building connectivities: " << endl;
  printf("Tranposition CellPoint -> PointCell\n");
  auto begin = std::chrono::steady_clock::now();

  con_PointCell = mesh::transpose( con_CellPoint );

  auto end = std::chrono::steady_clock::now();
  std::cout << "Time elapsed= = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl << std::endl;;
  begin = std::chrono::steady_clock::now();
  printf("Intersection CellPoint, PointCell -> CellCell\n");

  con_CellCell = mesh::intersect( con_CellPoint, con_PointCell );

  end = std::chrono::steady_clock::now();
  std::cout << "Time elapsed= = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl << std::endl;
  begin = std::chrono::steady_clock::now();
  printf("Build (Facet, Points), (Cells, facet)\n");

  tie(con_FacetPoint, con_CellFacet) = mesh::build( refCellEl.dim-1, con_CellPoint, con_CellCell );

  end = std::chrono::steady_clock::now();
  std::cout << "Time elapsed= = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl << std::endl;;

  con_FacetCell = mesh::transpose(con_CellFacet);

  // Build AABBs
  setAABBs();

}
