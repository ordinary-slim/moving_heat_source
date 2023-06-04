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
    adimR_tstep = 3
    adimR_domain = 4

    # read input
    problemInput = readInput( inputFile )

    fixedProblemInput = dict( problemInput )
    movingProblemInput = dict( problemInput )

    # Mesh
    meshInputFixed, meshInputMoving = {}, {}
    meshInputFixed["points"], meshInputFixed["cells"], meshInputFixed["cell_type"] = mesh(boxDomain)
    meshInputMoving["points"], meshInputMoving["cells"], meshInputMoving["cell_type"] = meshAroundHS(adimR_domain, movingProblemInput)

    meshFixed  = mhs.Mesh(meshInputFixed)
    meshMoving = mhs.Mesh(meshInputMoving)

    # Problem params
    # set dt
    dt = setAdimR( adimR_tstep, fixedProblemInput )
    for input in [fixedProblemInput, movingProblemInput,]:
        input["dt"] = dt

    #set MRF business NO TRANSPORT
    movingProblemInput["isAdvection"] = 1
    movingProblemInput["advectionSpeedX"] = -fixedProblemInput["HeatSourceSpeedX"]
    movingProblemInput["speedFRF_X"]      = fixedProblemInput["HeatSourceSpeedX"]
    movingProblemInput["HeatSourceSpeedX"] = 0.0

    pFixed         = Problem(meshFixed, fixedProblemInput, caseName="fixed")
    pFRF           = Problem(meshFixed, fixedProblemInput, caseName="FRF")
    pMoving        = Problem(meshMoving, movingProblemInput, caseName="moving")

    maxIter = pFixed.input["maxIter"]

    for iteration in range(maxIter):
        pFRF.iterate()

        pFixed.setAssembling2External( True )
        pMoving.setAssembling2External( True )

        # PRE-ITERATE AND DOMAIN OPERATIONS
        pMoving.domain.resetActivation()
        pFixed.domain.resetActivation()
        pMoving.intersectExternal( pFixed, False, False )
        pFixed.preiterate(False)
        pMoving.preiterate(False)
        pMoving.intersectExternal( pFixed, False, True )
        pFixed.substractExternal( pMoving, False, True )

        #Dirichet gamma
        pFixed.setGamma2Dirichlet()

        # Pre-assembly, updating free dofs
        pMoving.preAssemble(True)
        pFixed.preAssemble(True)
        # Allocate linear system
        ls = mhs.LinearSystem( pMoving, pFixed )
        ls.cleanup()

        pMoving.assemble()
        pFixed.assemble()

        pFixed.assembleDirichletGamma( pMoving )
        pMoving.assembleNeumannGamma( pFixed )

        ls.assemble()

        ls.solve()

        pFixed.gather()
        pMoving.gather()

        # Get inactive points information from other
        pFixed.unknown.interpolateInactive( pMoving.unknown )
        pMoving.unknown.interpolateInactive( pFixed.unknown )

        pFixed.postIterate()
        pMoving.postIterate()

        pFRF.writepos(
                )
        pFixed.writepos(
            nodeMeshTags={ "gammaNodes":pFixed.gammaNodes, },
                )
        pMoving.writepos(
            nodeMeshTags={ "gammaNodes":pMoving.gammaNodes, },
            )
