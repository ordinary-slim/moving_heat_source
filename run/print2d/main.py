import MovingHeatSource as mhs
import numpy as np
import meshzoo
import sys
import pdb

def mesh(box, meshDen=1, variant="zigzag", cell_type="triangle3"):
    '''
    Variant = "zigzag",  or "up", "down", "center"
    '''
    nelsX = int(meshDen*(box[1]-box[0]))
    nelsY = int(meshDen*(box[3]-box[2]))
    if cell_type=="triangle3":
        points, cells = meshzoo.rectangle_tri(
            np.linspace(box[0], box[1], nelsX+1),
            np.linspace(box[2], box[3], nelsY+1),
            #cell_type=cell_type,
            variant=variant)
    elif cell_type=="quad4":
        points, cells = meshzoo.rectangle_quad(
            np.linspace(box[0], box[1], nelsX+1),
            np.linspace(box[2], box[3], nelsY+1),
            cell_type=cell_type
            #variant="zigzag",  # or "up", "down", "center"
        )
    else:
        exit()
    cells = cells.astype( np.uint32 )
    return points, cells, cell_type

def meshAroundHS( adimR, problemInput, meshDen=4, cell_type="triangle3" ):
    radius = problemInput["radius"]
    initialPositionX = problemInput["initialPositionX"]
    initialPositionY = problemInput["initialPositionY"]
    trailLength = adimR * radius
    capotLength = min( trailLength, 3*radius )
    halfLengthY = min( trailLength, capotLength )
    box = [initialPositionX - trailLength, initialPositionX + capotLength,
           initialPositionY - halfLengthY, initialPositionY + halfLengthY,
           ]

    return mesh(box, meshDen=meshDen, cell_type=cell_type)

def deactivateBelowSurface(p, surfaceZ = 0):
    nels = p.domain.mesh.nels
    activeEls = []
    for ielem in range( nels ):
        e = p.domain.mesh.getElement( ielem )
        if (e.getCentroid()[1] < surfaceZ):
            activeEls.append( ielem )
    substrateEls = mhs.MeshTag( p.domain.mesh, p.domain.mesh.dim, activeEls )
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
    try:
        runReference = (sys.argv[1]=="True")
    except IndexError:
        pass
    inputFile = "input.yaml"
    boxDomain = [-25, 25, -5, 1]
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
    meshInputFixed["points"], meshInputFixed["cells"], meshInputFixed["cell_type"] = mesh(boxDomain, meshDen=meshDen, cell_type="quad4")
    ##meshDen = 4
    meshInputMoving["points"], meshInputMoving["cells"], meshInputMoving["cell_type"] = meshAroundHS(adimR_domain, movingProblemInput, meshDen=meshDen, cell_type="quad4")

    meshFixed  = mhs.Mesh(meshInputFixed)
    meshMoving = mhs.Mesh(meshInputMoving)

    # mhs.Problem params
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

    pFixed         = mhs.Problem(meshFixed, fixedProblemInput, caseName="fixed")
    pFRF           = mhs.Problem(meshFixed, fixedProblemInput, caseName="FRF")
    pFineFRF       = None
    if runReference:
        pFineFRF = mhs.Problem(meshFixed, referenceProblemInput, caseName="FineFRF")
    pMoving        = mhs.Problem(meshMoving, movingProblemInput, caseName="moving")

    Tfinal = pFixed.input["Tfinal"]

    tol = 1e-7

    for p in [pFixed, pFRF]:
        deactivateBelowSurface( p )
    if pFineFRF:
        deactivateBelowSurface( pFineFRF )


    # Set up printer
    mdwidth = 0.99
    mdheight = 1.99
    printerFRF = mhs.Printer( pFRF, mdwidth, mdheight )

    while (pFRF.time < Tfinal - tol) :
        # Setup print
        p1 = pFRF.mhs.position
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
            p1 = pFineFRF.mhs.position
            p2 = p1 + pFineFRF.mhs.speed * finedt
            # Print
            printerFineFRF.deposit( p1, p2, pFineFRF.domain.activeElements )

            # FRF ITERATE
            pFineFRF.iterate()

            pFineFRF.writepos()

    printerMoving = mhs.Printer( pMoving, mdwidth, mdheight )
    activeElsFixed = mhs.MeshTag( pFixed.domain.activeElements )

    while (pFixed.time < Tfinal - tol) :
        p1 = np.array(pFixed.mhs.position)
        # Put this in a loop
        # MY SCHEME ITERATE
        # PRE-ITERATE AND DOMAIN OPERATIONS
        pMoving.domain.resetActivation()
        pFixed.domain.setActivation(activeElsFixed)

        pMoving.intersectExternal(pFixed, updateGamma=False)

        pFixed.preiterate( canPreassemble=False )
        pMoving.preiterate( canPreassemble=False )

        pMoving.intersectExternal(pFixed, updateGamma=False)
        pFixed.substractExternal(pMoving, updateGamma=True)
        pMoving.updateInterface( pFixed )
        # Print
        # Setup print
        p2 = np.array(pFixed.mhs.position)
        p2 = p2 - 0.001*(p2 - p1)
        printerMoving.deposit( p1, p2, pMoving.domain.activeElements )

        #Dirichet gamma
        pFixed.setGamma2Dirichlet()
        # Pre-assembly, updating free dofs
        pMoving.preAssemble(allocateLs=True)
        pFixed.preAssemble(allocateLs=True)
        ls = mhs.LinearSystem.Create( self.pMoving, self.pFixed )
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
        pFixed.uniteExternal(pMoving, updateGamma=False)
        activeElsFixed = mhs.MeshTag( pFixed.domain.activeElements )
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
