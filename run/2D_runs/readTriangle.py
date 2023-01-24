import sys
sys.path.insert(1, '..')
sys.path.insert(1, '../../Debug/')
import MovingHeatSource as mhs
from readInput import *
import numpy as np
import meshzoo
import pdb

# two triangles and one quad
points, cells = meshzoo.rectangle_quad(
    np.linspace(0.0, 10.0, 3),
    np.linspace(0.0, 10.0, 3),
    cell_type="quad4",  # or "quad8", "quad9"
)

d = formatInputFile( "input.txt" )
d = parseInput(d)
d["points"] = points
d["cells"] = cells
d["cell_type"]="quad4"

print( "points:", points )
print( "tpoints:", type(points) )
print( "cells:", cells )
print( "tcells:", type(cells) )

p = mhs.Problem()
p.initialize( d )
print( p.solution )
