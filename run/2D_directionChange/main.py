import MovingHeatSource as mhs
from MovingHeatSource.gcode import gcode2laserPath
import numpy as np
import meshzoo
from AdaptiveStepper import AdaptiveStepper
import pdb

Tfinal = 2
inputFile = "input.yaml"
tol = 1e-7
adimR_domain = 10 
problemInput = mhs.readInput( inputFile )
adimR = problemInput["radius"] / np.linalg.norm(problemInput["HeatSourceSpeed"])

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

def meshSquare(center, sideLength, meshDen=4):
    cell_type="quad4"
    box = [center[0] - sideLength / 2,
           center[0] + sideLength / 2,
           center[1] - sideLength / 2,
           center[1] + sideLength / 2]
    return meshBox(box, meshDen=meshDen)

def meshAroundHS( adimR, problemInput, meshDen=4 ):
    radius = problemInput["radius"]
    initialPosition = problemInput["initialPosition"]
    trailLength = adimR * radius
    capotLength = min( trailLength, 2*radius )
    halfLengthY = min( trailLength, capotLength )
    box = [initialPosition[0] - trailLength, initialPosition[0] + capotLength,
           initialPosition[1] - halfLengthY, initialPosition[1] + halfLengthY]
    return meshBox(box, meshDen)

if __name__=="__main__":
    boxPhys = [-10, 10, -10, 10]

    # read input
    fixedProblemInput = dict( problemInput )
    movingProblemInput = dict( problemInput )

    FRFInput = dict( problemInput )

    # Mesh
    meshInputPhys, meshInputMoving = {}, {}
    meshInputPhys["points"], meshInputPhys["cells"], meshInputPhys["cell_type"] = meshBox(boxPhys)
    meshInputMoving["points"], meshInputMoving["cells"], meshInputMoving["cell_type"] = meshAroundHS(adimR_domain, movingProblemInput)
    #meshInputMoving["points"], meshInputMoving["cells"], meshInputMoving["cell_type"] = meshSquare(movingProblemInput["initialPosition"], 2*adimR_domain, meshDen=8)

    meshFixed           = mhs.Mesh(meshInputPhys)
    meshMoving           = mhs.Mesh(meshInputMoving)

    movingProblemInput["isAdvection"] = 1
    movingProblemInput["advectionSpeed"] = -fixedProblemInput["HeatSourceSpeed"]
    movingProblemInput["speedFRF"]      = fixedProblemInput["HeatSourceSpeed"]
    movingProblemInput["HeatSourceSpeed"] = np.zeros(3)

    pFRF     = mhs.Problem(meshFixed, FRFInput, caseName="FRF")
    pFixed   = mhs.Problem(meshFixed, fixedProblemInput, caseName="fixed")
    pMoving  = mhs.Problem(meshMoving, movingProblemInput, caseName="moving")


    maxIter = pFRF.input["maxIter"]

    pFRF.setDt( 1*adimR )
    pFixed.setDt( 1*adimR )
    pMoving.setDt( 1*adimR )

    # Set path
    for p in [pFRF, pFixed]:
        p.mhs.setPath( *gcode2laserPath( "Path.gcode" ) )

    myDriver = AdaptiveStepper( pFixed, pMoving,
                               adimMaxSubdomainSize=adimR_domain, rotateSubdomain=True )

    ## RIGHT
    '''
    while not(pFRF.mhs.path.isOver):
        pFRF.iterate()
        pFRF.writepos()
    '''

    while not(pFixed.mhs.path.isOver(pFixed.time)):
        myDriver.iterate()
