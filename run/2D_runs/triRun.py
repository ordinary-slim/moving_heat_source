import sys
sys.path.insert(1, '..')
sys.path.insert(1, '../../Debug/')
import MovingHeatSource as mhs
from readInput import *
import numpy as np
import meshzoo
import meshio

points, cells = meshzoo.rectangle_tri(
    np.linspace(-5, 5, 11),
    np.linspace(-5, 5, 11),
    variant="zigzag",  # or "up", "down", "center"
)

d = formatInputFile( "input.txt" )
d = parseInput(d)
d["points"] = points
d["cells"] = cells
d["cell_type"]="triangle3"

p = mhs.Problem()
p.initialize( d )

for iteration in range( d["maxIter"] ):
    p.iterate()


#export
tmpCells = [
    ("triangle", cells),
        ]
mesh = meshio.Mesh(
    points,
    tmpCells,
    # Optionally provide extra data on points, cells, etc.
    point_data={"T": p.solution},
)

mesh.write(
    "foo.vtk",  # str, os.PathLike, or buffer/open file
    # file_format="vtk",  # optional if first argument is a path; inferred from extension
)
