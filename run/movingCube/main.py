import sys
sys.path.insert(1, '..')
sys.path.insert(1, '../../Release/')
import MovingHeatSource as mhs
import numpy as np
import meshzoo
from wrapper import Problem, readInput
import pdb

def mesh(box, meshDen=4):
    cell_type="hexa8"
    nelsX = int((box[1] - box[0])*meshDen)
    nelsY = int((box[3] - box[2])*meshDen)
    nelsZ = int((box[5] - box[4])*meshDen)

    points, cells = meshzoo.cube_hexa(
        np.linspace( box[0], box[1], nelsX+1),
        np.linspace( box[3], box[2], nelsY+1),
        np.linspace( box[5], box[4], nelsZ+1),
    )
    cells = cells.astype( np.uint32 )
    return points, cells, cell_type

def meshAroundHS( adimR, problemInput, meshDen=4 ):
    radius = problemInput["radius"]
    initialPositionX = problemInput["initialPositionX"]
    initialPositionY = problemInput["initialPositionY"]
    initialPositionZ = problemInput["initialPositionY"]
    trailLength = adimR * radius
    capotLength = min( trailLength, 3*radius )
    halfLengthY = min( trailLength, capotLength )
    halfLengthZ = halfLengthY
    box = [initialPositionX - trailLength, initialPositionX + capotLength,
           initialPositionY - halfLengthY, initialPositionY + halfLengthY,
           initialPositionZ - halfLengthZ, initialPositionZ + halfLengthZ,
           ]

    return mesh(box, meshDen)

def setAdimR( adimR, input ):
    r = input["radius"]
    HeatSourceSpeedX = max( abs(input["HeatSourceSpeedX"]), abs(input["advectionSpeedX"]))
    HeatSourceSpeedY = max( abs(input["HeatSourceSpeedY"]), abs(input["advectionSpeedY"]))
    HeatSourceSpeedZ = max( abs(input["HeatSourceSpeedZ"]), abs(input["advectionSpeedZ"]))
    speed  = np.linalg.norm( np.array( [HeatSourceSpeedX, HeatSourceSpeedY, HeatSourceSpeedZ] ) )
    return (adimR * r / speed)

if __name__=="__main__":
    runReference = False
    if (len(sys.argv)>1):
        if sys.argv[1]=="True":
            runReference = True
    inputFile = "input.txt"
    boxDomain = [-25, 25, -5, 5, -5, 0]
    adimR_tstep = 2
    adimR_domain = 4

    # read input
    problemInput = readInput( inputFile )

    fixedProblemInput = dict( problemInput )
    referenceProblemInput = dict( problemInput )
    movingProblemInput = dict( problemInput )

    # Mesh
    meshDen = 2
    meshInputFixed, meshInputMoving = {}, {}
    meshInputFixed["points"], meshInputFixed["cells"], meshInputFixed["cell_type"] = mesh(boxDomain, meshDen=meshDen)
    meshDen = 3
    meshInputMoving["points"], meshInputMoving["cells"], meshInputMoving["cell_type"] = meshAroundHS(adimR_domain, movingProblemInput, meshDen=meshDen)

    meshFixed  = mhs.Mesh(meshInputFixed)
    meshMoving = mhs.Mesh(meshInputMoving)

    # Problem params
    # set dt
    dt = setAdimR( adimR_tstep, fixedProblemInput )
    for input in [fixedProblemInput, movingProblemInput,]:
        input["dt"] = dt
    referenceProblemInput["dt"] = setAdimR( 0.2, referenceProblemInput )

    #set MRF business NO TRANSPORT
    movingProblemInput["isAdvection"] = 1
    movingProblemInput["advectionSpeedX"] = -fixedProblemInput["HeatSourceSpeedX"]
    movingProblemInput["speedFRF_X"]      = fixedProblemInput["HeatSourceSpeedX"]
    movingProblemInput["HeatSourceSpeedX"] = 0.0

    pFixed         = Problem(meshFixed, fixedProblemInput, caseName="fixed")
    pFRF           = Problem(meshFixed, fixedProblemInput, caseName="FRF")
    pFineFRF       = None
    if runReference:
        pFineFRF = Problem(meshFixed, referenceProblemInput, caseName="FineFRF")
    pMoving        = Problem(meshMoving, movingProblemInput, caseName="moving")

    Tfinal = pFixed.input["Tfinal"]

    tol = 1e-7
    while (pMoving.time < Tfinal - tol) :
        # FRF ITERATE
        pFRF.iterate()

        # Fine FRF ITERATE
        if pFineFRF:
            while (pFineFRF.time < pFRF.time - tol):
                pFineFRF.iterate()

        # MY SCHEME ITERATE
        pFixed.setAssembling2External( True )
        pMoving.setAssembling2External( True )
        # PRE-ITERATE AND DOMAIN OPERATIONS
        pMoving.domain.resetActivation()
        pFixed.domain.resetActivation()
        pMoving.intersectExternal( pFixed, False, False )
        pFixed.preiterate(False)
        pMoving.preiterate(False)
        pMoving.intersectExternal( pFixed, False, False )
        pFixed.substractExternal( pMoving, False, True )
        pMoving.updateInterface( pFixed )
        #Dirichet gamma
        pFixed.setGamma2Dirichlet()
        # Pre-assembly, updating free dofs
        pMoving.preAssemble(True)
        pFixed.preAssemble(True)
        # Allocate linear system
        ls = mhs.LinearSystem( pMoving, pFixed )
        ls.cleanup()
        # Assembly
        pMoving.assemble()
        pFixed.assemble()
        # Assembly Gamma
        pFixed.assembleDirichletGamma( pMoving )
        pMoving.assembleNeumannGamma( pFixed )
        # Build ls
        ls.assemble()
        # Solve ls
        ls.solve()
        # Recover solution
        pFixed.gather()
        pMoving.gather()
        pFixed.unknown.interpolateInactive( pMoving.unknown, False )
        pMoving.unknown.interpolateInactive( pFixed.unknown, True )
        # Post iteration
        pFixed.postIterate()
        pMoving.postIterate()

        # POSTPROCESSING ALL
        pFRF.writepos(
                )
        if pFineFRF:
            pFineFRF.writepos(
                    )

        pFixed.writepos(
            nodeMeshTags={ "gammaNodes":pFixed.gammaNodes, },
                )
        pMoving.writepos(
            nodeMeshTags={ "gammaNodes":pMoving.gammaNodes, },
            )