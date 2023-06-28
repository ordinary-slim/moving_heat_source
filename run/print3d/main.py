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

def deactivateBelowSurface(p, surfaceZ = 0):
    nels = p.domain.mesh.nels
    activeEls = []
    for ielem in range( nels ):
        e = p.domain.mesh.getElement( ielem )
        if (e.getCentroid()[2] < surfaceZ):
            activeEls.append( ielem )
    substrateEls = mhs.mark( p.domain.mesh, p.domain.mesh.dim, activeEls )
    p.domain.setActivation( substrateEls )

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
    boxDomain = [-25, 25, -5, 5, -5, 1]
    adimR_tstep = 0.5
    adimR_domain = 2

    # read input
    problemInput = readInput( inputFile )

    fixedProblemInput = dict( problemInput )
    referenceProblemInput = dict( problemInput )
    movingProblemInput = dict( problemInput )

    # Mesh
    meshDen = 2
    meshInputFixed, meshInputMoving = {}, {}
    meshInputFixed["points"], meshInputFixed["cells"], meshInputFixed["cell_type"] = mesh(boxDomain, meshDen=meshDen)
    meshDen = 2
    meshInputMoving["points"], meshInputMoving["cells"], meshInputMoving["cell_type"] = meshAroundHS(adimR_domain, movingProblemInput, meshDen=meshDen)

    meshFixed  = mhs.Mesh(meshInputFixed)
    meshMoving = mhs.Mesh(meshInputMoving)

    # Problem params
    # set dt
    dt = setAdimR( adimR_tstep, fixedProblemInput )
    for input in [fixedProblemInput, movingProblemInput,]:
        input["dt"] = dt
    finedt = setAdimR( 0.5, referenceProblemInput )
    referenceProblemInput["dt"] = finedt

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

    for p in [pFixed, pFRF]:
        deactivateBelowSurface( p )
    if pFineFRF:
        deactivateBelowSurface( pFineFRF )


    # Set up printer
    mdwidth = 1
    mdheight = 2
    printerFRF = mhs.Printer( pFRF, mdwidth, mdheight )

    '''
    while (pFRF.time < Tfinal - tol) :
        # Setup print
        p1 = pFRF.mhs.currentPosition
        p2 = p1 + pFRF.mhs.speed * dt
        # Print
        printerFRF.deposit( p1, p2, pFRF.domain.activeElements )

        # FRF ITERATE
        pFRF.iterate()

        pFRF.writepos(
                )

    if runReference:
        printerFineFRF = mhs.Printer( pFineFRF, mdwidth, mdheight )

        while (pFineFRF.time < Tfinal - tol) :
            # Setup print
            p1 = pFineFRF.mhs.currentPosition
            p2 = p1 + pFineFRF.mhs.speed * finedt
            # Print
            printerFineFRF.deposit( p1, p2, pFineFRF.domain.activeElements )

            # FRF ITERATE
            pFineFRF.iterate()

            pFineFRF.writepos()
    '''

    printerMoving = mhs.Printer( pMoving, mdwidth, mdheight )
    activeElsFixed = mhs.MeshTag( pFixed.domain.activeElements )

    while (pFixed.time < Tfinal - tol) :
        p1 = np.array(pFixed.mhs.currentPosition)
        # Put this in a loop
        # MY SCHEME ITERATE
        pFixed.setAssembling2External( True )
        pMoving.setAssembling2External( True )
        # PRE-ITERATE AND DOMAIN OPERATIONS
        pMoving.domain.resetActivation()
        pFixed.domain.setActivation(activeElsFixed)

        pMoving.intersectExternal( pFixed, False )

        pFixed.preiterate(False)
        pMoving.preiterate(False)

        pMoving.intersectExternal( pFixed, False )
        pFixed.substractExternal( pMoving, True )
        pMoving.updateInterface( pFixed )
        # Print
        # Setup print
        p2 = np.array(pFixed.mhs.currentPosition)
        printerMoving.deposit( p1, p2, pMoving.domain.activeElements )

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

        # Union
        pFixed.uniteExternal( pMoving, False )
        activeElsFixed = mhs.MeshTag( pFixed.domain.activeElements )
        pMoving.unknown.interpolateInactive( pFixed.unknown, True )

        #DEBUGGING
        prevValFixed = mhs.Function(pFixed.previousValues[0])
        prevValMoving = mhs.Function(pMoving.previousValues[0])
        #DEBUGGING
        # Post iteration
        pFixed.postIterate()
        pMoving.postIterate()

        pFixed.writepos(
            functions={"prevVal":prevValFixed},
            nodeMeshTags={ "gammaNodes":pFixed.gammaNodes, },
            cellMeshTags={ "elsOwnedByOther":pFixed.elsOwnedByOther, },
                )
        pMoving.writepos(
            functions={"prevVal":prevValMoving},
            nodeMeshTags={ "gammaNodes":pMoving.gammaNodes, },
            cellMeshTags={ "elsOwnedByOther":pMoving.elsOwnedByOther, },
            )
