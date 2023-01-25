import sys
sys.path.insert(1, '..')
sys.path.insert(1, '../../Debug/')
import MovingHeatSource as mhs
from readInput import *
import numpy as np
import meshzoo
import meshio

cell_type="quad4"
points, cells = meshzoo.rectangle_quad(
    np.linspace(-5, 5, 11),
    np.linspace(-5, 5, 11),
    cell_type=cell_type
    #variant="zigzag",  # or "up", "down", "center"
)

d = formatInputFile( "input.txt" )
d = parseInput(d)
d["points"] = points
d["cells"] = cells
d["cell_type"]=cell_type

print( "Points = \n", points )
p = mhs.Problem()
p.initialize( d )

for iteration in range( d["maxIter"] ):
    p.iterate()

el = p.mesh.getElement(0)
K = np.zeros((el.nnodes,el.nnodes))
for inode in range(el.nnodes):
    for jnode in range(el.nnodes):
        for igp in range(el.ngpoints):
            K[inode, jnode] += np.dot( el.GradBaseGpVals[inode][igp], el.GradBaseGpVals[jnode][igp] ) * el.vol *el.gpweight[igp]

#export
tmpCells = [
    ("quad", cells),
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
