import MovingHeatSource as mhs
import numpy as np
import meshzoo
import pdb

nsteps=5

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

def run():
    inputFile = "input.py"
    boxDomain = [-0.005, 0.005, -0.025, +0.025, 0.0, 0.02]

    # read input
    problemInput = mhs.readInput( inputFile )
    mdwidth = 2.5e-3
    mdheight = 2.5e-3
    problemInput["heatSourceHeight"] = mdheight
    problemInput["heatSourceWidth"] = mdwidth

    fixedProblemInput = dict( problemInput )

    # Mesh
    meshDen = 100
    meshInputFixed = {}
    meshInputFixed["points"], meshInputFixed["cells"], meshInputFixed["cell_type"] = mesh(boxDomain, meshDen=meshDen)

    meshFixed  = mhs.Mesh(meshInputFixed)

    # mhs.Problem params
    # set dt
    dt = 1
    fixedProblemInput["dt"] = dt

    p           = mhs.Problem(meshFixed, fixedProblemInput, caseName="case")

    for p in [p]:
        deactivateBelowSurface( p, surfaceZ=0.01 )

    # Set up printer
    printerFRF = mhs.Printer( p, mdwidth, mdheight )

    for _ in range(nsteps):
        # Setup print
        p1 = p.mhs.currentPosition
        p2 = p1 + p.mhs.speed * dt
        # Print
        printerFRF.deposit( p1, p2, p.domain.activeElements )
        p.mhs.markHeatedElements(p1, p2)
        # Apply convection to all external boundaries
        p.setConvection()

        p.iterate()

        p.writepos()

def test():
    run()
    # COMPARISON
    allFilesSame = True
    for ifile in range(nsteps):
        newds = "post_case/case_{}.vtu".format( ifile+1 )
        refds = "post_case_reference/case_{}.vtu".format( ifile+1 )
        if not(mhs.meshio_comparison( newds, refds )):
               allFilesSame = False
               break
    assert allFilesSame

if __name__=="__main__":
    test()
