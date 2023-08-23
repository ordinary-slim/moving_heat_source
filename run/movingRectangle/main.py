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
    initialPosition = problemInput["initialPosition"]
    trailLength = adimR * radius
    capotLength = min( trailLength, 3*radius )
    halfLengthY = min( trailLength, capotLength )
    box = [initialPosition[0] - trailLength, initialPosition[0] + capotLength,
           initialPosition[1] - halfLengthY, initialPosition[1] + halfLengthY]
    return mesh(box, meshDen)

def setAdimR( adimR, input ):
    r = input["radius"]
    speed = max( np.linalg.norm( input["HeatSourceSpeed"] ), np.linalg.norm( input["advectionSpeed"] ) )
    return (adimR * r / speed)

if __name__=="__main__":
    runReference = False
    if (len(sys.argv)>1):
        if sys.argv[1]=="True":
            runReference = True
    inputFile = "input.yaml"
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

    movingProblemInput["isAdvection"] = 1
    movingProblemInput["advectionSpeed"] = -fixedProblemInput["HeatSourceSpeed"]
    movingProblemInput["speedFRF"]      = fixedProblemInput["HeatSourceSpeed"]
    movingProblemInput["HeatSourceSpeed"] = np.zeros(3)

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
        # PRE-ITERATE AND DOMAIN OPERATIONS
        pMoving.domain.resetActivation()
        pFixed.domain.resetActivation()
        pMoving.intersectExternal(pFixed, updateGamma=False)
        pFixed.preiterate( canPreassemble=False )
        pMoving.preiterate( canPreassemble=False )
        pMoving.intersectExternal(pFixed, updateGamma=False)
        pFixed.substractExternal(pMoving, updateGamma=True)
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
