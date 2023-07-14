import MovingHeatSource as mhs
import numpy as np
import meshzoo
import sys

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
    trailLength = adimR * radius
    capotLength = min( trailLength, 3*radius )
    halfLengthY = min( trailLength, capotLength )
    box = [initialPositionX - trailLength, initialPositionX + capotLength,
           initialPositionY - halfLengthY, initialPositionY + halfLengthY]
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
    inputFile = "input.py"
    boxDomain = [-25, 25, -5, 5]
    adimR_tstep = 10
    adimR_domain = 11

    # read input
    problemInput = mhs.readInput( inputFile )

    fixedProblemInput = dict( problemInput )
    referenceProblemInput = dict( problemInput )
    movingProblemInput = dict( problemInput )

    # Mesh
    meshDen = 2
    meshInputFixed, meshInputMoving = {}, {}
    meshInputFixed["points"], meshInputFixed["cells"], meshInputFixed["cell_type"] = mesh(boxDomain, meshDen=meshDen)
    meshInputMoving["points"], meshInputMoving["cells"], meshInputMoving["cell_type"] = meshAroundHS(adimR_domain, movingProblemInput, meshDen=meshDen)

    meshFixed  = mhs.Mesh(meshInputFixed)
    meshMoving = mhs.Mesh(meshInputMoving)

    # mhs.Problem params
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

    pFixed         = mhs.Problem(meshFixed, fixedProblemInput, caseName="fixed")
    pFRF           = mhs.Problem(meshFixed, fixedProblemInput, caseName="FRF")
    if runReference:
        pFineFRF = mhs.Problem(meshFixed, referenceProblemInput, caseName="FineFRF")
    pMoving        = mhs.Problem(meshMoving, movingProblemInput, caseName="moving")

    Tfinal = pFixed.input["Tfinal"]

    tol = 1e-7
    while (pMoving.time < Tfinal - tol) :
        # FRF ITERATE
        pFRF.iterate()

        # Fine FRF ITERATE
        if runReference:
            while (pFineFRF.time < pFRF.time - tol):
                pFineFRF.iterate()

        # MY SCHEME ITERATE
        pFixed.setAssembling2External( True )
        pMoving.setAssembling2External( True )
        # PRE-ITERATE AND DOMAIN OPERATIONS
        pMoving.domain.resetActivation()
        pFixed.domain.resetActivation()
        pMoving.intersectExternal( pFixed, False )
        pFixed.preiterate(False)
        pMoving.preiterate(False)
        pMoving.intersectExternal( pFixed, False )
        pFixed.substractExternal( pMoving, False)
        pMoving.updateInterface( pFixed )
        #Dirichet gamma
        pFixed.setGamma2Dirichlet()
        # Pre-assembly, updating free dofs
        pMoving.preAssemble(True)
        pFixed.preAssemble(True)
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
        if runReference:
            pFineFRF.writepos(
                    )

        pFixed.writepos(
            nodeMeshTags={ "gammaNodes":pFixed.gammaNodes,
                          },
                )
        pMoving.writepos(
            nodeMeshTags={ "gammaNodes":pMoving.gammaNodes,
                          },
            )
