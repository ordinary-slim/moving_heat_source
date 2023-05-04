import sys
sys.path.insert(1, '..')
sys.path.insert(1, '../../Debug/')
import MovingHeatSource as mhs
import numpy as np
import meshzoo
from wrapper import Problem, readInput
import pdb

maxIter = 40
dt = 0.5
Tfinal = 20
initialT = 10

def mesh(box):
    cell_type="quad4"
    meshDen = 1
    points, cells = meshzoo.rectangle_quad(
        np.linspace(box[0], box[1], meshDen*(box[1]-box[0])+1),
        np.linspace(box[2], box[3], meshDen*(box[3]-box[2])+1),
        cell_type=cell_type
        #variant="zigzag",  # or "up", "down", "center"
    )
    cells = cells.astype( int )
    return points, cells, cell_type

if __name__=="__main__":
    inputFile = "input.txt"
    box = [-16, 16, -5, 5]

    # Read input
    problemInput = readInput( inputFile )

    # Mesh
    meshInput = {}
    meshInput["points"], meshInput["cells"], meshInput["cell_type"] = mesh(box)
    m = mhs.Mesh(meshInput)

    # Initialize problems
    p  = Problem(m, problemInput, caseName="case")


    f = lambda pos : initialT
    p.forceState( f )

    it = 0
    time = 0.0
    while (( time < Tfinal) and ( it < maxIter )):
        time += dt
        it += 1
        p.iterate()
        p.writepos()
