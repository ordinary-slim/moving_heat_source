import MovingHeatSource as mhs
from MovingHeatSource.gcode import gcode2laserPath
import numpy as np
import meshzoo
from AdaptiveStepper import AdaptiveStepper
import pdb

Tfinal = 2
inputFile = "input.yaml"
tol = 1e-7
adimR_domain = 3
problemInput = mhs.readInput( inputFile )
adimR = problemInput["radius"] / problemInput["HeatSourceSpeedX"]

def meshBox(box, meshDen=4):
    cell_type="quad4"
    nelsX = int(meshDen*(box[1]-box[0])) +1
    nelsY = int(meshDen*(box[3]-box[2])) +1
    points, cells = meshzoo.rectangle_quad(
        np.linspace(box[0], box[1], nelsX),
        np.linspace(box[2], box[3], nelsY),
        cell_type=cell_type
        #variant="zigzag",  # or "up", "down", "center"
    )
    cells = cells.astype( int )
    return points, cells, cell_type

def iterate(p):
    p.iterate()

if __name__=="__main__":
    boxPhys = [-10, 10, -10, 10]

    # read input
    FRFInput = dict( problemInput )

    # Mesh
    meshInputPhys, meshInputMoving = {}, {}
    meshInputPhys["points"], meshInputPhys["cells"], meshInputPhys["cell_type"] = meshBox(boxPhys, meshDen=2)

    meshFixed           = mhs.Mesh(meshInputPhys)

    pFRF     = mhs.Problem(meshFixed, FRFInput, caseName="FRF")

    maxIter = pFRF.input["maxIter"]

    pFRF.setDt( 1*adimR )

    # Set path
    pFRF.mhs.setPath( *gcode2laserPath( "Path.gcode" ) )

    while not(pFRF.mhs.path.isOver):
        iterate( pFRF )
