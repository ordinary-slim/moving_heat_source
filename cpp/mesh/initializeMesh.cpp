#include <chrono>
#include "Mesh.h"
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

namespace py = pybind11;

mesh::Mesh::Mesh(const py::dict &input) {
  //SET ELEMENT TYPE
  ElementType cell_type_flag =
    string2ElementType( py::cast<string>( input["cell_type"] ) );

  int ngpointsCell = -1, ngpointsFacet = -1;
  if (input.contains("numberOfGaussPoints")){
    ngpointsCell = py::cast<int>( input["numberOfGaussPointsCells"] );
  }
  if (input.contains("numberOfGaussPointsFacets")){
    ngpointsFacet = py::cast<int>( input["numberOfGaussPointsFacets"] );
  }
  // reference element. no support for mixed meshes yet
  refCellEl = ReferenceElement(cell_type_flag, ngpointsCell);
  ElementType FacetElType = getFacetElType(cell_type_flag);
  refFacetEl = ReferenceElement(FacetElType, ngpointsFacet);

  //READ POINTS
  py::array points = input["points"];
  nnodes = points.shape(0);
  if (input.contains("dimension")){
    dim = py::cast<int>( input["dimension"] );
  } else {
    dim    = points.shape(1);
  }
  pos.resize( nnodes, 3 );
  pos.setZero();
  auto aux_points = points.unchecked<double>();
  for ( int ipoint = 0; ipoint < nnodes; ipoint++) {
    for ( int idim = 0; idim < dim; idim++) {
      pos(ipoint, idim) =  aux_points(ipoint, idim);
    }
  }
  //READ ELEMENTS
  py::array cells    = input["cells"];
  nels          = cells.shape(0);
  nnodes_per_el = cells.shape(1);
  con_CellPoint.con.resize( nels );
  elementTypes.resize( nels );
  auto aux_cells = cells.unchecked<unsigned int>();
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

  tie(con_FacetPoint, con_CellFacet) = mesh::buildBoundaryConnectivities( con_CellPoint, con_CellCell );

  end = std::chrono::steady_clock::now();
  std::cout << "Time elapsed= = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl << std::endl;;

  con_FacetCell = mesh::transpose(con_CellFacet);

  // Build AABBs
  buildAABBTree();

}
