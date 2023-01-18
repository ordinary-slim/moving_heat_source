import meshzoo
import meshio
import numpy as np
import pdb

points, cells = meshzoo.rectangle_quad(
    np.linspace(0.0, 10.0, 11),
    np.linspace(0.0, 10.0, 11),
    cell_type="quad4",  # or "quad8", "quad9"
)
tmpCells = [
    ("quad", cells),
        ]

# compute fake field
T = np.zeros( len( points ) )
for idx, p in enumerate(points):
    T[idx] = p[0]**2 + p[1]**2

mesh = meshio.Mesh(
    points,
    tmpCells,
    # Optionally provide extra data on points, cells, etc.
    point_data={"T": T},
)

mesh.write(
    "foo.vtk",  # str, os.PathLike, or buffer/open file
    # file_format="vtk",  # optional if first argument is a path; inferred from extension
)

#meshio.write_points_cells( "foo.vtk", points, cells,)
