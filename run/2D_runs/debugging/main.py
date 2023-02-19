import sys
sys.path.insert(1, '..')
sys.path.insert(1, '../../../Debug/')
import MovingHeatSource as mhs
import numpy as np
import meshzoo
from wrapper import Problem
import math
import pdb

def mesh(box, meshDen=1):
    cell_type="quad4"
    points, cells = meshzoo.rectangle_quad(
        np.linspace(box[0], box[1], math.ceil(meshDen*(box[1]-box[0])+1)),
        np.linspace(box[2], box[3], math.ceil(meshDen*(box[3]-box[2])+1)),
        cell_type=cell_type
        #variant="zigzag",  # or "up", "down", "center"
    )
    cells = cells.astype( int )
    return points, cells, cell_type

if __name__=="__main__":
    inputFile = "input.txt"
    bgBox = [-10, 10, -5, +5]
    adimR = 1

    bgProblem = Problem("background")
    fgProblem = Problem("foreground")
    # INITIALIZATIONS
    # Mesh BG
    points, cells, cell_type = mesh(bgBox)
    bgProblem.setMesh( points, cells, cell_type )

    for p in [fgProblem, bgProblem]:
        p.parseInput( inputFile )

    # Mesh FG
    center = np.array( [bgProblem.input["initialPositionX"],
                        bgProblem.input["initialPositionY"],
                        bgProblem.input["initialPositionZ"]], dtype=np.float64,)
    lSide  = 2.5*bgProblem.input["radius"]
    fgBox  = [center[0] - lSide/2,
              center[0] + lSide/2,
              center[1] - lSide/2,
              center[1] + lSide/2,]
    points, cells, cell_type = mesh(fgBox, meshDen=4)
    fgProblem.setMesh( points, cells, cell_type )

    # Set advection FG
    fgProblem.input["isAdvection"] = 1
    fgProblem.input["advectionSpeedX"] = -bgProblem.input["speedX"]
    fgProblem.input["speedX"] = 0
    fgProblem.input["power"] = 0

    for p in [fgProblem, bgProblem]:
        p.initialize()

    # ITERATE BG ONCE
    bgProblem.updateFRFpos()#get tn+1 positions (not tn)
    bgProblem.iterate()#assembly + solve
    bgProblem.writepos()

    fgProblem.getFromExternal( bgProblem.mesh, bgProblem.unknown, bgProblem.shiftFRF )
    fgProblem.writepos()

    # MARCH FG, INTERPOLATE BACK TO BG
    for it in range(fgProblem.input["maxIter"]):
        fgProblem.updateFRFpos()
        #MARCH
        fgProblem.postIterate()
        bgProblem.postIterate()

        bgProblem.getFromExternal( fgProblem.mesh, fgProblem.unknown, fgProblem.shiftFRF )

        fgProblem.writepos()
        bgProblem.writepos()
