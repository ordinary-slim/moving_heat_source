import sys
sys.path.insert(1, '..')
sys.path.insert(1, '../../Release/')
import MovingHeatSource as mhs
import numpy as np
import meshzoo
from wrapper import Problem, readInput
import pdb

def mesh(box, meshDen=4):
    cell_type="quad4"
    nelsX = int( meshDen*(box[1]-box[0])+1)
    nelsY = int( meshDen*(box[3]-box[2])+1)
    points, cells = meshzoo.rectangle_quad(
        np.linspace(box[0], box[1], nelsX),
        np.linspace(box[2], box[3], nelsY),
        cell_type=cell_type
        #variant="zigzag",  # or "up", "down", "center"
    )
    cells = cells.astype( int )
    return points, cells, cell_type

def meshAroundHS( adimR, problemInput, meshDen=4 ):
    radius = problemInput["radius"]
    initialPositionX = problemInput["initialPositionX"]
    initialPositionY = problemInput["initialPositionY"]
    halfLength = adimR * radius
    box = [initialPositionX - halfLength, initialPositionX + halfLength,
           initialPositionY - halfLength, initialPositionY + halfLength]
    return mesh(box, meshDen)

def setAdimR( adimR, input ):
    r = input["radius"]
    HeatSourceSpeedX = max( abs(input["HeatSourceSpeedX"]), abs(input["advectionSpeedX"]))
    HeatSourceSpeedY = max( abs(input["HeatSourceSpeedY"]), abs(input["advectionSpeedY"]))
    HeatSourceSpeedZ = max( abs(input["HeatSourceSpeedZ"]), abs(input["advectionSpeedZ"]))
    speed  = np.linalg.norm( np.array( [HeatSourceSpeedX, HeatSourceSpeedY, HeatSourceSpeedZ] ) )
    return (adimR * r / speed)

if __name__=="__main__":
    inputFile = "input.txt"
    boxDomain = [-16, 16, -5, 5]
    adimR = 1

    # read input
    problemInput = readInput( inputFile )

    fixedProblemInput = dict( problemInput )
    movingProblemInput = dict( problemInput )

    # Mesh
    meshInputFixed, meshInputMoving = {}, {}
    meshInputFixed["points"], meshInputFixed["cells"], meshInputFixed["cell_type"] = mesh(boxDomain)
    meshInputMoving["points"], meshInputMoving["cells"], meshInputMoving["cell_type"] = meshAroundHS(adimR, movingProblemInput)

    meshFixed  = mhs.Mesh(meshInputFixed)
    meshMoving = mhs.Mesh(meshInputMoving)

    # Problem params
    # set dt
    dt = setAdimR( adimR, fixedProblemInput )
    for input in [fixedProblemInput, movingProblemInput,]:
        input["dt"] = dt

    #set MRF business NO TRANSPORT
    movingProblemInput["isAdvection"] = 1
    movingProblemInput["advectionSpeedX"] = -fixedProblemInput["HeatSourceSpeedX"]
    movingProblemInput["speedFRF_X"]      = fixedProblemInput["HeatSourceSpeedX"]
    movingProblemInput["HeatSourceSpeedX"] = 0.0

    pFixed         = Problem(meshFixed, fixedProblemInput, caseName="fixed")
    pMoving        = Problem(meshMoving, movingProblemInput, caseName="moving")

    maxIter = pFixed.input["maxIter"]
    # FORWARD
    for iteration in range(maxIter):

        pFixed.setAssembling2External( True )
        pMoving.setAssembling2External( True )

        # PRE-ITERATE
        pFixed.preiterate()
        pMoving.preiterate()
        pFixed.updateFRFpos()
        pMoving.updateFRFpos()

        # DOMAIN OPERATIONS
        pFixed.substractExternal( pMoving, True, True )

        #pFixed.gather()
        #pMoving.gather()

        #pFixed.postIterate()
        #pMoving.postIterate()

        pFixed.writepos(
            nodeMeshTags={ "gammaNodes":pFixed.gammaNodes, },
                )
        pMoving.writepos()
