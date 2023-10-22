import MovingHeatSource as mhs
import numpy as np
import meshzoo
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
    initialPosition = problemInput["initialPosition"]
    trailLength = adimR * radius
    capotLength = min( trailLength, 3*radius )
    halfLengthY = min( trailLength, capotLength )
    halfLengthZ = halfLengthY
    box = [initialPosition[0] - trailLength, initialPosition[0] + capotLength,
           initialPosition[1] - halfLengthY, initialPosition[1] + halfLengthY,
           initialPosition[2] - halfLengthZ, initialPosition[2] + halfLengthZ,
           ]

    return mesh(box, meshDen)

def setAdimR( adimR, input ):
    r = input["radius"]
    speed = max( np.linalg.norm( input["HeatSourceSpeed"] ), np.linalg.norm( input["advectionSpeed"] ) )
    return (adimR * r / speed)

if __name__=="__main__":
    runReference = False
            runReference = True
    inputFile = "input.yaml"
    boxDomain = [-25, 25, -5, 5, -5, 0]
    adimR_tstep = 2
    adimR_domain = 4

    # read input
    problemInput = mhs.readInput( inputFile )

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

    # mhs.Problem params
    # set dt
    dt = setAdimR( adimR_tstep, fixedProblemInput )
    for input in [fixedProblemInput, movingProblemInput,]:
        input["dt"] = dt
    referenceProblemInput["dt"] = setAdimR( 0.2, referenceProblemInput )

    #set MRF business NO TRANSPORT
    movingProblemInput["isAdvection"] = 1
    movingProblemInput["advectionSpeed"] = -fixedProblemInput["HeatSourceSpeed"]
    movingProblemInput["speedFRF"]      = fixedProblemInput["HeatSourceSpeed"]
    movingProblemInput["HeatSourceSpeed"] = np.zeros(3)

    pFixed         = mhs.Problem(meshFixed, fixedProblemInput, caseName="fixed")
    pFRF           = mhs.Problem(meshFixed, fixedProblemInput, caseName="FRF")
    pFineFRF       = None
    if runReference:
        pFineFRF = mhs.Problem(meshFixed, referenceProblemInput, caseName="FineFRF")
    pMoving        = mhs.Problem(meshMoving, movingProblemInput, caseName="moving")

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
        # PRE-ITERATE AND DOMAIN OPERATIONS
        pMoving.domain.resetActivation()
        pFixed.domain.resetActivation()
        pMoving.intersectExternal( pFixed, False, False )
        pFixed.preIterate( canPreassemble=False )
        pMoving.preIterate( canPreassemble=False )
        pMoving.intersectExternal( pFixed, False, False )
        pFixed.substractExternal( pMoving, False, True )
        pMoving.updateInterface( pFixed )
        #Dirichet gamma
        pFixed.setGamma2Dirichlet()
        pMoving.setGamma2Neumann()
        # Pre-assembly, updating free dofs
        pMoving.preAssemble(allocateLs=True)
        pFixed.preAssemble(allocateLs=True)
        ls = mhs.LinearSystem.Create( self.pMoving, self.pFixed )
        # Assembly
        pMoving.assemble( pFixed )
        pFixed.assemble( pMoving )
        # Build ls
        ls.assemble()
        # Solve ls
        ls.solve()
        # Recover solution
        pFixed.gather()
        pMoving.gather()
        pFixed.unknown.interpolateInactive( pMoving.unknown, ignoreOutside = False )
        pMoving.unknown.interpolateInactive( pFixed.unknown, ignoreOutside = True )
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
