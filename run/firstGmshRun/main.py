import MovingHeatSource as mhs
from MovingHeatSource.gcode import gcode2laserPath
import numpy as np
import meshio
import meshzoo
from AdaptiveStepper import AdaptiveStepper

inputFile = "input.yaml"
problemInput = mhs.readInput( inputFile )
tol = 1e-7

def deactivateBelowSurface(p,
                           surfacey = 0):
    nels = p.domain.mesh.nels
    activeels = []
    for ielem in range( nels ):
        e = p.domain.mesh.getElement( ielem )
        if (e.getCentroid()[1] < surfacey):
            activeels.append( ielem )
    substrateels = mhs.MeshTag( p.domain.mesh, p.domain.mesh.dim, activeels )
    p.domain.setActivation( substrateels )

def readMesh(gmshFile):
    m = meshio.read( gmshFile )
    mDict = {}
    mDict["points"] = m.points

    cells = np.vstack( [cell.data for cell in m.cells] )
    mDict["cells"] = cells
    mDict["cell_type"] = "quad4"

    return mhs.Mesh( mDict )

if __name__=="__main__":
    # read input
    fixedProblemInput = dict( problemInput )

    # Mesh
    meshFixed = readMesh("untitled.msh")

    pFixed   = mhs.Problem(meshFixed, fixedProblemInput, caseName="fixed")
    pFRF   = mhs.Problem(meshFixed, fixedProblemInput, caseName="FRF")

    maxIter = pFixed.input["maxIter"]


    # Set path, deactivate
    for p in [pFRF, pFixed]:
        p.mhs.setPath( *gcode2laserPath( "Path.gcode" ) )
        deactivateBelowSurface( p )

    myDriver = AdaptiveStepper( pFixed,
                               adimMaxSubdomainSize=10,
                               rotateSubdomain=True )


    # Set up printer
    y_width = 0.99
    myDriver.setupPrinter( y_width, 1 )

    while not(pFixed.mhs.path.isOver(pFixed.time)):
        myDriver.iterate()


    printerFRF = mhs.Printer( pFRF, y_width, 1 )
    while not(pFRF.mhs.path.isOver(pFRF.time)):
        # Setup print
        p1 = pFRF.mhs.position
        p2 = pFRF.mhs.path.interpolatePosition( pFRF.time + pFRF.dt )
        # Print
        printerFRF.deposit( p1, p2, pFRF.domain.activeElements )
        pFRF.iterate()
        pFRF.writepos()

