import MovingHeatSource as mhs
import numpy as np
import meshzoo
import sys

inputFile = "input.yaml"
boxDomain = [-25, 25, -5, 5, -5, 1]
problemInput = mhs.readInput( inputFile )
Tfinal = problemInput["Tfinal"]
tol = 1e-7
mdwidth = 1
mdheight = 1

# read input
problemInput = mhs.readInput( inputFile )


def mesh(box, elementSize=0.25):
    cell_type="hexa8"
    nelsX = int((box[1] - box[0])/elementSize)
    nelsY = int((box[3] - box[2])/elementSize)
    nelsZ = int((box[5] - box[4])/elementSize)

    points, cells = meshzoo.cube_hexa(
        np.linspace( box[0], box[1], nelsX+1),
        np.linspace( box[3], box[2], nelsY+1),
        np.linspace( box[5], box[4], nelsZ+1),
    )
    cells = cells.astype( np.uint32 )
    return points, cells, cell_type

def meshAroundHS( adimR, problemInput, elementSize=0.25 ):
    radius = problemInput["radius"]
    initialPosition = problemInput["initialPosition"]
    trailLength = adimR * radius
    capotLength = min( trailLength, 2*radius )
    halfLengthY = min( trailLength, capotLength )
    halfLengthZ = halfLengthY
    box = [initialPosition[0] - trailLength, initialPosition[0] + capotLength,
           initialPosition[1] - halfLengthY, initialPosition[1] + halfLengthY,
           initialPosition[2] - halfLengthZ, initialPosition[2] + halfLengthZ,
           ]

    meshDict = {}
    meshDict["points"], meshDict["cells"], meshDict["cell_type"] = mesh(box, elementSize)
    meshDict["dimension"] = 3
    return mhs.Mesh( meshDict )

def deactivateBelowSurface(p, surfaceZ = 0):
    nels = p.domain.mesh.nels
    activeEls = []
    for ielem in range( nels ):
        e = p.domain.mesh.getElement( ielem )
        if (e.getCentroid()[2] < surfaceZ):
            activeEls.append( ielem )
    substrateEls = mhs.MeshTag( p.domain.mesh, p.domain.mesh.dim, activeEls )
    p.domain.setActivation( substrateEls )

def setAdimR( adimR, input ):
    r = input["radius"]
    speed  = np.linalg.norm( input["HeatSourceSpeed"] )
    return (adimR * r / speed)

def getMesh( boxDomain, elementSize ):
    meshInputFixed, meshInputMoving = {}, {}
    meshInputFixed["points"], meshInputFixed["cells"], meshInputFixed["cell_type"] = mesh(boxDomain, elementSize=elementSize)
    return mhs.Mesh(meshInputFixed)

def runReference():
    elementSize = 0.5
    mesh = getMesh( boxDomain, elementSize )
    pFRF           = mhs.Problem(mesh, dict(problemInput), caseName="FRF")

    pFRF.setDt( setAdimR(0.5, pFRF.input ) )

    deactivateBelowSurface( pFRF )

    printerFRF = mhs.Printer( pFRF, mdwidth, mdheight/2, mdheight/2 )

    while (pFRF.time < Tfinal - tol) :
        # Setup print
        p1 = pFRF.mhs.position
        p2 = p1 + pFRF.mhs.speed * pFRF.dt
        # Print
        printerFRF.deposit( p1, p2, pFRF.domain.activeElements )

        # FRF ITERATE
        pFRF.iterate()

        pFRF.writepos(
                )

def runCoupled():
    adimR_tstep = 2
    adimR_domain = 4
    fixedProblemInput = dict( problemInput )
    movingProblemInput = dict( problemInput )

    # Mesh
    elementSize = 0.5
    meshFixed = getMesh( boxDomain, elementSize )
    meshMoving = meshAroundHS(adimR_domain, movingProblemInput, elementSize=elementSize)

    # mhs.Problem params
    # set dt
    dt = setAdimR( adimR_tstep, fixedProblemInput )
    for input in [fixedProblemInput, movingProblemInput]:
        input["dt"] = dt

    #set MRF business NO TRANSPORT
    movingProblemInput["isAdvection"] = 1
    movingProblemInput["advectionSpeed"] = -fixedProblemInput["HeatSourceSpeed"]
    movingProblemInput["speedDomain"]      = fixedProblemInput["HeatSourceSpeed"]
    movingProblemInput["HeatSourceSpeed"] = np.zeros(3)

    pFixed         = mhs.Problem(meshFixed, fixedProblemInput, caseName="fixed")
    pMoving         = mhs.Problem(meshMoving, movingProblemInput, caseName="moving")

    deactivateBelowSurface( pFixed )

    # Set up printer

    printerMoving = mhs.Printer( pMoving, mdwidth, mdheight/2, mdheight/2 )
    physicalEls = mhs.MeshTag( pFixed.domain.activeElements )

    while (pFixed.time < Tfinal - tol) :
        p1 = np.array(pFixed.mhs.position)
        # Put this in a loop
        # MY SCHEME ITERATE
        # PRE-ITERATE AND DOMAIN OPERATIONS
        pMoving.domain.resetActivation()
        pFixed.domain.setActivation(physicalEls)

        pMoving.intersectExternal(pFixed, updateGamma=False)

        pFixed.preiterate( canPreassemble=False )
        pMoving.preiterate( canPreassemble=False )

        pMoving.intersectExternal(pFixed, updateGamma=False)
        pFixed.substractExternal(pMoving, updateGamma=True)
        pMoving.updateInterface( pFixed )
        # Print
        # Setup print
        p2 = np.array(pFixed.mhs.position)
        printerMoving.deposit( p1, p2, pMoving.domain.activeElements )

        #Dirichet gamma
        pFixed.setGamma2Dirichlet()
        pMoving.setGamma2Neumann()
        # Pre-assembly, updating free dofs
        pMoving.preAssemble(allocateLs=True)
        pFixed.preAssemble(allocateLs=True)
        ls = mhs.LinearSystem.Create( pMoving, pFixed )
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

        # Union
        pFixed.uniteExternal(pMoving, updateGamma=False)
        physicalEls = mhs.MeshTag( pFixed.domain.activeElements )
        pMoving.unknown.interpolateInactive( pFixed.unknown, ignoreOutside = True )

        # Post iteration
        pFixed.postIterate()
        pMoving.postIterate()

        pFixed.writepos(
            nodeMeshTags={ "gammaNodes":pFixed.gammaNodes, },
                )
        pMoving.writepos(
            nodeMeshTags={ "gammaNodes":pMoving.gammaNodes, },
            )

if __name__=="__main__":
    isRunReference = ("--run-reference" in sys.argv)
    isOnlyRunReference = ("--only-reference"  in sys.argv)
    if isRunReference or isOnlyRunReference:
        runReference()
    if not(isOnlyRunReference):
        runCoupled()
