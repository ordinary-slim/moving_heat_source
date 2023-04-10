import sys
sys.path.insert(1, '..')
sys.path.insert(1, '../../Release/')
import MovingHeatSource as mhs
import numpy as np
import meshzoo
from wrapper import Problem
import pdb

maxIter = 50
dt = 0.1
Tfinal = 5
initialT = 10

def mesh(box):
    cell_type="quad4"
    meshDen = 2
    points, cells = meshzoo.rectangle_quad(
        np.linspace(box[0], box[1], meshDen*(box[1]-box[0])+1),
        np.linspace(box[2], box[3], meshDen*(box[3]-box[2])+1),
        cell_type=cell_type
        #variant="zigzag",  # or "up", "down", "center"
    )
    cells = cells.astype( int )
    return points, cells

if __name__=="__main__":
    inputFile = "input.txt"
    box = [-16, 16, -5, 5]

    p             = Problem("case")

    points, cells = mesh(box)
    for p in [p,]:
        p.input["points"] = points
        p.input["cells"] = cells
        p.input["cell_type"]="quad4"

    #read input
    for p in [p,]:
        p.parseInput( inputFile )

    for p in [p,]:
        p.initialize()

    f = lambda pos : initialT
    p.forceState( f )

    it = 0
    time = 0.0
    while (( time < Tfinal) and ( it < maxIter )):
        time += dt
        it += 1
        p.iterate()
        p.writepos()
