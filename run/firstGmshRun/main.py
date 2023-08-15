import MovingHeatSource as mhs
from MovingHeatSource.gcode import gcode2laserPath
import numpy as np
import meshio
import meshzoo
from AdaptiveStepper import AdaptiveStepper
import pdb

inputFile = "input.yaml"
tol = 1e-7
adimR_domain = 10 
problemInput = mhs.readInput( inputFile )
adimR = problemInput["radius"] / np.linalg.norm(problemInput["HeatSourceSpeed"])

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


def meshBox(box, meshDen=4):
    cell_type="quad4"
    nelsX = int(meshDen*(box[1]-box[0])) +1
    nelsY = int(meshDen*(box[3]-box[2])) +1
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
    capotLength = min( trailLength, 2*radius )
    halfLengthY = min( trailLength, capotLength )
    box = [initialPosition[0] - trailLength, initialPosition[0] + capotLength,
           initialPosition[1] - halfLengthY, initialPosition[1] + halfLengthY]
    return meshBox(box, meshDen)

if __name__=="__main__":
    # read input
    fixedProblemInput = dict( problemInput )
    movingProblemInput = dict( problemInput )

    # Mesh
    meshInputPhys, meshInputMoving = {}, {}
    meshInputMoving["points"], meshInputMoving["cells"], meshInputMoving["cell_type"] = meshAroundHS(adimR_domain, movingProblemInput)

    meshFixed = readMesh("untitled.msh")
    meshMoving = mhs.Mesh(meshInputMoving)

    movingProblemInput["isAdvection"] = 1
    movingProblemInput["advectionSpeed"] = -fixedProblemInput["HeatSourceSpeed"]
    movingProblemInput["speedFRF"]      = fixedProblemInput["HeatSourceSpeed"]
    movingProblemInput["HeatSourceSpeed"] = np.zeros(3)

    pFixed   = mhs.Problem(meshFixed, fixedProblemInput, caseName="fixed")
    pMoving  = mhs.Problem(meshMoving, movingProblemInput, caseName="moving")

    maxIter = pFixed.input["maxIter"]


    # Set path, deactivate
    for p in [pFixed]:
        p.mhs.setPath( *gcode2laserPath( "Path.gcode" ) )
        deactivateBelowSurface( p )

    myDriver = AdaptiveStepper( pFixed,
                               pMoving,
                               adimMaxSubdomainSize=adimR_domain,
                               rotateSubdomain=True )


    # Set up printer
    mdwidth = 1
    mdheight = 1
    myDriver.setupPrinter( mdwidth, mdheight )

    while not(pFixed.mhs.path.isOver(pFixed.time)):
        myDriver.iterate()
